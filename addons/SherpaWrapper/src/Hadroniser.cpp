/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2026  Laurent Forthomme
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <CepGen/Core/Exception.h>
#include <CepGen/Core/RunParameters.h>
#include <CepGen/Event/Event.h>
#include <CepGen/Modules/EventModifierFactory.h>
#include <CepGen/Modules/RandomGeneratorFactory.h>
#include <CepGen/Physics/Hadroniser.h>
#include <CepGen/Physics/PDG.h>
#include <CepGen/Utils/RandomGenerator.h>
#include <CepGen/Utils/String.h>

#include "CepGenSherpa/DiffractiveState.h"
#include "CepGenSherpa/SherpaFramework.h"

using namespace std::string_literals;

namespace cepgen::sherpa {
  /// Interface to the Sherpa fragmentation (AHADIC++ cluster model by default) and hadron decays (HADRONS++)
  /// algorithms for the dissociated beam systems
  class Hadroniser final : public hadr::Hadroniser {
  public:
    explicit Hadroniser(const ParametersList& params)
        : hadr::Hadroniser(params),
          random_generator_(RandomGeneratorFactory::get().build(steer<ParametersList>("randomGenerator"))),
          primordial_kt_(steer<double>("primordialKt")),
          tolerance_(steer<double>("momentumTolerance")),
          fragmentation_(steer<std::string>("fragmentation")),
          hadron_decays_(steer<bool>("hadronDecays")),
          run_card_(steer<std::string>("runCard")),
          verbosity_(steer<int>("verbosity")) {
      if (utils::toLower(fragmentation_) == "none")
        throw CG_FATAL("sherpa:Hadroniser") << "A fragmentation algorithm is required.";
    }

    static ParametersDescription description() {
      auto desc = hadr::Hadroniser::description();
      desc.setDescription("Interface to the Sherpa fragmentation and hadron decays of dissociated beam systems");
      desc.add("maxTrials", 10).setDescription("maximum number of hadronisation attempts per dissociated system");
      desc.add("randomGenerator", RandomGeneratorFactory::get().describeParameters("stl"))
          .setDescription("random number generator to use for the partonic content of the dissociated systems");
      desc.add("primordialKt", 0.)
          .setDescription("RMS of the spectator transverse momentum relative to the beam axis (GeV)");
      desc.add("momentumTolerance", 1.e-6)
          .setDescription("relative tolerance on the four-momentum balance between the system and its products");
      desc.add("fragmentation", "Ahadic"s)
          .setDescription("Sherpa fragmentation algorithm (e.g. 'Ahadic' cluster model, or 'Pythia8' Lund string)");
      desc.add("hadronDecays", true).setDescription("decay the unstable hadrons with the Sherpa HADRONS++ module?");
      desc.add("runCard", ""s).setDescription("optional Sherpa YAML run card for the fragmentation/decays tuning");
      desc.add("verbosity", 0).setDescription("Sherpa output level (0 = errors only, 2 = information)");
      return desc;
    }

    /// Parse a Sherpa setting
    void readString(const std::string& parameter) override {
      if (sherpa_initialised_)
        throw CG_FATAL("sherpa:Hadroniser") << "Cannot apply '" << parameter << "': Sherpa is already set up.";
      const auto setting = utils::trim(parameter);
      const auto key_end = setting.find_first_of(":=");
      if (key_end == std::string::npos || key_end == 0)
        throw CG_FATAL("sherpa:Hadroniser")
            << "Invalid parameter '" << parameter << "'. Expecting 'KEY: value' or 'KEY=value'.";
      const auto key = utils::toUpper(utils::trim(setting.substr(0, key_end)));
      for (const auto& protected_key :
           {"EVENT_TYPE", "EVENTS", "FRAGMENTATION", "RANDOM_SEED", "RUNDATA", "INIT_ONLY", "PATH", "SHERPA_CPP_PATH"})
        if (key == protected_key)
          throw CG_FATAL("sherpa:Hadroniser")
              << "Sherpa setting '" << key << "' is steered by the CepGen interface, and cannot be overridden.";
      if (key == "HADRON_DECAYS" && utils::toUpper(setting).find("MODEL") != std::string::npos)
        throw CG_FATAL("sherpa:Hadroniser") << "Hadron decays model is steered through the 'hadronDecays' parameter.";
      settings_.emplace_back(setting);
    }

    void initialise() override {
      mass_threshold_ = PDG::get().mass(PDG::proton) + PDG::get().mass(PDG::piZero);  // lightest two-body final state
      roles_.clear();
      const auto mode = runParameters().kinematics().incomingBeams().mode();
      if (mode == mode::Kinematics::InelasticElastic || mode == mode::Kinematics::InelasticInelastic)
        roles_[Particle::Role::IncomingBeam1] = Particle::Role::OutgoingBeam1;
      if (mode == mode::Kinematics::ElasticInelastic || mode == mode::Kinematics::InelasticInelastic)
        roles_[Particle::Role::IncomingBeam2] = Particle::Role::OutgoingBeam2;
    }

    bool run(Event& event, double& weight, bool fast) override {
      weight = 1.;
      if (!fragment_remnants_)
        return true;
      // also vetoed during the integration, for these points not to contribute to the cross section
      for (const auto& [incoming, outgoing] : roles_)
        if (event.oneWithRole(outgoing).momentum().mass() < mass_threshold_)
          return false;  // do not start hadronisation if diffractive mass under the threshold
      if (fast)
        return true;

      auto& framework = Framework::get();
      if (!sherpa_initialised_) {  // see configuration()
        framework.initialise(configuration());
        sherpa_initialised_ = true;
      }

      std::unordered_map<Particle::Role, hadr::Fragments> hadronised;  // nothing is written until all systems succeeded
      int num_trials = 0;
      {
        std::lock_guard<std::mutex> lock(framework.mutex());  // Sherpa algorithms rely on global objects
        for (const auto& [incoming, outgoing] : roles_) {
          if (event(outgoing).size() != 1)  // already fragmented by another algorithm
            continue;
          const auto& beam = event.oneWithRole(incoming);
          const auto diffractive_seed = sherpa::DiffractiveState(
              *random_generator_, beam, event.oneWithRole(outgoing).momentum(), primordial_kt_);
          if (!diffractive_seed.valid()) {
            CG_WARNING("sherpa:Hadroniser") << "Invalid remnants kinematics for " << outgoing << " (mass "
                                            << event.oneWithRole(outgoing).momentum().mass() << " GeV).";
            return false;
          }
          hadr::Fragments fragments;
          int trials_for_system = 0;
          do {
            if (trials_for_system++ >= max_trials_)
              return false;
            ++num_trials;
          } while (!diffractive_seed.hadronise(fragments, tolerance_));
          hadronised[outgoing] = std::move(fragments);
        }
      }
      event.metadata["sherpa:num_hadronisation_trials"] = num_trials;

      for (const auto& [role, fragments] : hadronised)
        update(event, role, fragments);
      return true;
    }

  private:
    /// Build Sherpa configuration from the steered parameters
    /// \note Only done for first hadronisation
    Framework::Configuration configuration() const {
      Framework::Configuration configuration;
      configuration.arguments = {"EVENT_TYPE: HadronDecay",  // no hard process/MPI
                                 "EVENTS: 1",
                                 "FRAGMENTATION: " + fragmentation_,
                                 "HADRON_DECAYS: {Model: "s + (hadron_decays_ ? "HADRONS++" : "None") + "}"};
      if (seed_ > 0)
        configuration.arguments.emplace_back("RANDOM_SEED: " + std::to_string(seed_));
      configuration.arguments.insert(configuration.arguments.end(), settings_.begin(), settings_.end());
      configuration.run_card = run_card_;
      std::ostringstream defaults;
      defaults << "BEAMS: [11, -11]\n"  // electron beams -> no PDF
               << "BEAM_ENERGIES: 50.\n"
               << "PDF_LIBRARY: None\n"
               << "OUTPUT: " << verbosity_ << "\n"
               << "WRITE_REFERENCES_FILE: false\n";
      configuration.defaults = defaults.str();  // user settings
      return configuration;
    }

    /// Add the hadrons to the event, once the full event was successfully hadronised
    static void update(Event& event, Particle::Role role, const hadr::Fragments& fragments) {
      const auto system_id = event.oneWithRole(role).id();
      const auto system_momentum = event[system_id].momentum();
      event[system_id].children().clear();

      std::unordered_map<int, int> identifiers;  // fragment index -> CepGen identifier
      for (const auto& fragment : fragments) {
        auto& particle = event.addParticle(role).get();
        particle.setIntegerPdgId(fragment.pdg_id)
            .setStatus(fragment.stable ? Particle::Status::FinalState : Particle::Status::Resonance)
            .setMomentum(Momentum(fragment.momentum).lorentzBoost(system_momentum), /*off_shell=*/true);
        const auto mother = identifiers.find(fragment.mother);
        particle.addMother(event[mother != identifiers.end() ? mother->second : system_id]);
        identifiers[fragment.index] = particle.id();
      }
      event[system_id].setStatus(Particle::Status::Fragmented);
    }

    const std::unique_ptr<utils::RandomGenerator> random_generator_;
    const double primordial_kt_;
    const double tolerance_;
    const std::string fragmentation_;
    const bool hadron_decays_;
    const std::string run_card_;
    const int verbosity_;

    std::unordered_map<Particle::Role, Particle::Role> roles_;  ///< (incoming, outgoing) dissociating beams
    std::vector<std::string> settings_;                         ///< User-steered Sherpa settings
    double mass_threshold_{0.};                                 ///< Lowest mass allowing a two-body decay
    bool sherpa_initialised_{false};
  };
}  // namespace cepgen::sherpa
using SherpaHadroniser = cepgen::sherpa::Hadroniser;
REGISTER_MODIFIER("sherpa", SherpaHadroniser);
