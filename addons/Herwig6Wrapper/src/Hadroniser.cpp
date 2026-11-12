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
#include <CepGen/Physics/Kinematics.h>
#include <CepGen/Physics/Modes.h>
#include <CepGen/Physics/PDG.h>
#include <CepGen/Utils/Message.h>
#include <CepGen/Utils/RandomGenerator.h>
#include <CepGen/Utils/StreamCollector.h>
#include <CepGen/Utils/String.h>

#include <cstring>

#include "CepGenHerwig6/DiffractiveState.h"
#include "CepGenHerwig6/Herwig6Interface.h"

namespace cepgen::herwig6 {
  namespace {  // static part of the algorithm, as Herwig 6 is Fortran, and dominated with common blocks
    std::unordered_map<std::string, double*> realParameters() {
      return {{"B1LIM", &hwpram_.b1lim},  // b meson-recombination parameter (default = 0)
              {"BTCLM", &hwpram_.btclm},  // multiplier to mass parameter to adjust splitting/grouping
                                          // of spectator partons into colour-singlets (default = 1)

              // fission threshold mass parameterisation (M_f = CLMAX**CLPOW + (m_q1 + m_q2)**CLPOW
              {"CLMAX", &hwpram_.clmax},         // maximum cluster mass parameter in cluster hadronisation/fission
                                                 // model (default = 3.35 GeV)
              {"CLPOW", &hwpram_.clpow},         // power-dependence cluster mass parameter
              {"CLSMR(1)", &hwpram_.clsmr[0]},   // cluster mass smearing parameter (no b quark)
              {"CLSMR(2)", &hwpram_.clsmr[1]},   // cluster mass smearing parameter (b quark)
              {"PSPLT(1)", &hwpram_.psplt[0]},   // steepness of mass spectrum (d, u, s quarks)
              {"PSPLT(2)", &hwpram_.psplt[1]}};  // steepness of mass spectrum (b quark)
    }
    std::unordered_map<std::string, int*> integerParameters() {
      return {{"CLDIR(1)", &hwpram_.cldir[0]},  // cluster dir. memory switch during non-pert. hadr. (d, u, s quarks)
              {"CLDIR(2)", &hwpram_.cldir[1]},  // cluster dir. memory switch during non-pert. hadr. (b quark)
              {"IOPREM", &hwpram_.ioprem}};  // Gauss. mass spectrum (0)/child as soft remnant cluster (1) (default = 1)
    }
    /// Parameter-dependent state of the (process-wide) Herwig 6 common blocks; empty until it is set up
    ParametersList& herwigConfiguration() {
      static ParametersList configuration;
      return configuration;
    }
  }  // namespace

  /// Interface to the Herwig 6 cluster hadronisation algorithm for the dissociated beam systems
  class Hadroniser final : public cepgen::hadr::Hadroniser {
  public:
    explicit Hadroniser(const ParametersList& params)
        : cepgen::hadr::Hadroniser(params),
          random_generator_(RandomGeneratorFactory::get().build(steer<ParametersList>("randomGenerator"))),
          momentum_tolerance_(steer<double>("momentumTolerance")),
          primordial_kt_(steer<double>("primordialKt")) {}

    // no call to the Herwig 6 library here, as this is also used for the documentation generation
    static ParametersDescription description() {
      auto desc = cepgen::hadr::Hadroniser::description();
      desc.setDescription("Interface to the Herwig 6 cluster hadronisation of dissociated beam systems");
      desc.add("maxTrials", 10).setDescription("maximum number of hadronisation attempts per dissociated system");
      desc.add("randomGenerator", RandomGeneratorFactory::get().describeParameters("stl"))
          .setDescription("random number generator to use for the partonic content of the dissociated systems");
      desc.add("momentumTolerance", 1.e-6)
          .setDescription("maximum fractional momentum tolerance to consider a hadronisation successful");
      desc.add("primordialKt", 0.)
          .setDescription("RMS of the spectator transverse momentum relative to the beam axis (GeV)");
      return desc;
    }

    void readString(const std::string& parameter) override {
      if (herwig_initialised_)
        throw CG_FATAL("herwig6:Hadroniser") << "Cannot apply '" << parameter << "': Herwig 6 is already set up.";
      const auto parts = utils::split(parameter, '=', true);
      if (parts.size() != 2)
        throw CG_FATAL("herwig6:Hadroniser") << "Invalid parameter '" << parameter << "'. Expecting 'NAME=value'.";
      const auto name = utils::toUpper(parts.at(0));
      if (realParameters().count(name) > 0) {
        settings_.set(name, std::stod(parts.at(1)));
      } else if (integerParameters().count(name) > 0)
        settings_.set(name, std::stoi(parts.at(1)));
      else
        throw CG_FATAL("herwig6:Hadroniser") << "Unsupported Herwig 6 parameter: '" << parts.at(0) << "'.";
    }

    void initialise() override {
      {  // collect the stdout output of Herwig 6 initialisation
        std::string output;
        utils::StreamCollector collector(output);
        herwig6::initialise();
      }
      if (hwprop_.idpdg[static_cast<int>(herwig6::DiffractiveState::InternalId::ud0)] != PDG::ud0 ||
          hwprop_.idpdg[static_cast<int>(herwig6::DiffractiveState::InternalId::proton)] != PDG::proton ||
          hwprop_.rmass[static_cast<int>(herwig6::DiffractiveState::InternalId::proton)] <= 0.)
        throw CG_FATAL("herwig6:Hadroniser")
            << "Herwig 6 particle table is not filled as expected.\n\t" << "Is the HWUDAT block data linked in?";
      if (!compatibleVersion())
        CG_WARNING("herwig6:Hadroniser") << "Common blocks layouts were not validated against linked Herwig 6 version.";
      threshold_ = hwprop_.rmass[static_cast<int>(herwig6::DiffractiveState::InternalId::proton)] +
                   hwprop_.rmass[static_cast<int>(herwig6::DiffractiveState::InternalId::piZero)];

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
        if (event.oneWithRole(outgoing).momentum().mass() < threshold_)
          return false;  // do not start hadronisation if diffractive mass under the threshold
      if (fast)
        return true;

      initialiseHerwig();

      std::unordered_map<Particle::Role, Fragments> hadronised;  // nothing is written until all systems succeeded
      int num_trials = 0;
      for (const auto& [incoming, outgoing] : roles_) {
        if (event(outgoing).size() != 1)  // already fragmented by another algorithm
          continue;
        const herwig6::DiffractiveState diffractive_seed(
            *random_generator_, event.oneWithRole(incoming), event.oneWithRole(outgoing).momentum(), primordial_kt_);
        auto& fragments = hadronised[outgoing];
        int trials_for_system = 0;
        do {
          if (trials_for_system++ >= max_trials_)
            return false;
          ++num_trials;
        } while (!diffractive_seed.hadronise(fragments, momentum_tolerance_));
      }
      event.metadata["herwig6:num_hadronisation_trials"] = num_trials;
      for (const auto& [role, fragments] : hadronised)
        update(event, role, fragments);
      return true;
    }

  private:
    /// Compute the Herwig 6 parameter-dependent constants
    /// \note This is only done when the first event is to be hadronised, as the post-initialisation configuration
    ///       blocks of the steering card are only parsed after the initialise() call
    void initialiseHerwig() {
      if (herwig_initialised_)  // already initialised
        return;

      constexpr double kReferenceBeamMomentum = 7000.;  // unused without a hard process, but has to be non-zero
      constexpr int kClusterOnlyProcess = 18000;  // no hard process, no shower tables, and no soft underlying event
      constexpr int kUnlimitedErrors = 2000000000;

      auto configuration = ParametersList{settings_}.set<unsigned long long>("seed", seed_);
      if (auto& current = herwigConfiguration(); !current.empty()) {
        if (current != configuration)
          throw CG_FATAL("herwig6:Hadroniser")
              << "Herwig 6 was already set up in this process with a different configuration.\n\t"
              << "Current: " << current << "\n\tRequested: " << configuration << ".";
        herwig_initialised_ = true;
        return;
      }
      std::memcpy(hwbmch_.part1, "P       ", 8);
      std::memcpy(hwbmch_.part2, "P       ", 8);
      hwproc_.pbeam1 = hwproc_.pbeam2 = kReferenceBeamMomentum;
      hwproc_.iproc = kClusterOnlyProcess;
      hwevnt_.maxer = kUnlimitedErrors;  // Herwig would otherwise stop the whole process after a few killed events
      hwevnt_.maxpr = 0;
      hwevnt_.nowgt = 0;  // do not search for the maximal weight of a hard process
      hwpram_.lrsud = hwpram_.lwsud = 0;
      hwpram_.iprint = 0;
      hwpram_.prvtx = 0;
      hwdspn_.syspin = 0;
      for (const auto& key : settings_.keys())  // map Herwig 6 common blocks to C++ configuration
        if (const auto reals = realParameters(); reals.count(key) > 0)
          settings_.fill(key, *reals.at(key));
        else
          settings_.fill(key, *integerParameters().at(key));
      if (seed_ > 0) {  // some numerology magic to split random number generator seed into 2 sub-parts
        hwevnt_.nrn[0] = static_cast<int>(seed_ % 2147483562) + 1;
        hwevnt_.nrn[1] = static_cast<int>((seed_ * 7) % 2147483398) + 1;
      }
      {  // start the actual Herwig 6 initialisation
        std::string output;
        utils::StreamCollector collector(output);
        hwuinc();
        hweini();
      }
      herwigConfiguration() = configuration;
      CG_DEBUG("herwig6:Hadroniser:initialiseHerwig").log([this](auto& log) {
        log << "Herwig 6 configuration: " << herwigConfiguration() << ". Fed blocks: dbl{";
        std::string sep;
        for (const auto& [key, value] : realParameters())
          log << sep << key << " -> " << *value, sep = ", ";
        sep.clear();
        log << ", int{";
        for (const auto& [key, value] : integerParameters())
          log << sep << key << " -> " << *value, sep = ", ";
        log << "}.";
      });
      herwig_initialised_ = true;
    }

    /// Add the hadrons to the event, once the full event was successfully hadronised
    void update(Event& event, Particle::Role role, const Fragments& fragments) const {
      // particles are appended to the vector the system lives in, hence only access it through its identifier
      const auto system_id = event.oneWithRole(role).id();
      const auto system_momentum = event[system_id].momentum();
      event[system_id].children().clear();  // events are recycled, with no reset of the system's filiation

      std::map<int, int> identifiers;  // event record position -> CepGen identifier
      for (const auto& fragment : fragments) {
        auto& particle = event.addParticle(role).get();
        particle.setIntegerPdgId(fragment.pdg_id)
            .setStatus(fragment.stable ? Particle::Status::FinalState : Particle::Status::Resonance)
            .setMomentum(Momentum(fragment.momentum).lorentzBoost(system_momentum), /*off_shell=*/true);
        const auto mother = identifiers.find(fragment.mother);
        particle.addMother(event[mother != identifiers.end() ? mother->second : system_id]);
        identifiers[fragment.index] = particle.id();
      }
      event[system_id].setStatus(Particle::Status::Fragmented);  // after the filiation, as it alters this status
    }

    const std::unique_ptr<utils::RandomGenerator> random_generator_;
    const double momentum_tolerance_;
    const double primordial_kt_;

    std::unordered_map<Particle::Role, Particle::Role> roles_;  ///< (incoming, outgoing) dissociating beams
    ParametersList settings_;                                   ///< User-steered Herwig 6 parameters
    double threshold_{0.};                                      ///< Lowest mass allowing a two-body decay
    bool herwig_initialised_{false};
  };
}  // namespace cepgen::herwig6
using Herwig6Hadroniser = cepgen::herwig6::Hadroniser;
REGISTER_MODIFIER("herwig6", Herwig6Hadroniser);
