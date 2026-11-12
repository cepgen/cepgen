/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2016-2026  Laurent Forthomme
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
#include <CepGen/Physics/Kinematics.h>
#include <CepGen/Physics/Modes.h>
#include <CepGen/Utils/RandomGenerator.h>
#include <Pythia8/Pythia.h>

#include <unordered_map>

#include "CepGenPythia8/DiffractiveState.h"
#include "CepGenPythia8/EventInterface.h"

using namespace std::string_literals;

namespace cepgen::pythia8 {
  /// Interface to the Pythia 8 hadronisation algorithm
  /// \note do not use full Pythia 8 event content -> but fragment beam remnants and decay central system afterwards.
  class Hadroniser final : public hadr::Hadroniser {
  public:
    explicit Hadroniser(const ParametersList& params)
        : hadr::Hadroniser(params),
          pythia_(new Pythia8::Pythia),
          random_generator_(RandomGeneratorFactory::get().build(steer<ParametersList>("randomGenerator"))),
          momentum_tolerance_(steer<double>("momentumTolerance")),
          primordial_kt_(steer<double>("primordialKt")),
          output_config_(steer<std::string>("outputConfig")) {}

    ~Hadroniser() override {
      if (!output_config_.empty())
        pythia_->settings.writeFile(output_config_, false);
    }

    static ParametersDescription description() {
      auto desc = hadr::Hadroniser::description();
      desc.setDescription("Interface to the Pythia 8 string hadronisation/fragmentation algorithm");
      desc.add("maxTrials", 10).setDescription("maximum number of hadronisation attempts per dissociated system");
      desc.add("randomGenerator", RandomGeneratorFactory::get().describeParameters("stl"))
          .setDescription("random number generator to use for the partonic content of the dissociated systems");
      desc.add("momentumTolerance", 1.e-6)
          .setDescription("maximum fractional momentum tolerance to consider a hadronisation successful");
      desc.add("primordialKt", 0.)
          .setDescription("RMS of the spectator transverse momentum relative to the beam axis (GeV)");
      desc.add("outputConfig", "last_pythia_config.cmd"s).setDescription("last Pythia configuration snapshot filename");
      return desc;
    }

    void readString(const std::string& param) override {
      if (!pythia_->readString(param))
        throw CG_FATAL("pythia8:Hadroniser") << "The Pythia 8 core failed to parse the following setting:\n\t" << param;
    }

    void initialise() override {
      // hadron-level standalone mode: no beams nor hard process, all systems are fed directly to the event record
      pythia_->settings.flag("ProcessLevel:all", false);
      if (seed_ == -1ll)
        pythia_->settings.flag("Random:setSeed", false);
      else {
        pythia_->settings.flag("Random:setSeed", true);
        pythia_->settings.mode("Random:seed", seed_);
      }
      if (!pythia_->init())
        throw CG_FATAL("pythia8:Hadroniser")
            << "Failed to initialise the Pythia 8 core!\n\t" << "See the message above for more details.";

      roles_.clear();
      const auto mode = runParameters().kinematics().incomingBeams().mode();
      if (mode == mode::Kinematics::InelasticElastic || mode == mode::Kinematics::InelasticInelastic)
        roles_[Particle::Role::IncomingBeam1] = Particle::Role::OutgoingBeam1;
      if (mode == mode::Kinematics::ElasticInelastic || mode == mode::Kinematics::InelasticInelastic)
        roles_[Particle::Role::IncomingBeam2] = Particle::Role::OutgoingBeam2;
    }

    bool run(Event& event, double& weight, bool fast) override {
      weight = 1.;

      // fragment the dissociated beams first, as nothing is written until all systems succeeded
      std::unordered_map<Particle::Role, hadr::Fragments> hadronised;
      if (!fast && fragment_remnants_) {  // do not treat beam remnants when running in fast mode
        int num_trials = 0;
        for (const auto& [incoming, outgoing] : roles_) {
          if (!event.hasRole(outgoing) || event(outgoing).size() != 1)  // already fragmented by another algorithm
            continue;
          const pythia8::DiffractiveState diffractive_state(*random_generator_,
                                                            event.oneWithRole(outgoing).momentum(),
                                                            event.oneWithRole(incoming),
                                                            *pythia_,
                                                            primordial_kt_);
          if (!diffractive_state.valid()) {
            CG_WARNING("pythia8:Hadroniser") << "Invalid remnants kinematics for " << outgoing << " (mass "
                                             << event.oneWithRole(outgoing).momentum().mass() << " GeV).";
            return false;
          }
          hadr::Fragments fragments;
          int trials_for_system = 0;
          do {
            if (trials_for_system++ >= max_trials_)
              return false;
            ++num_trials;
          } while (!diffractive_state.hadronise(fragments, momentum_tolerance_));
          hadronised[outgoing] = std::move(fragments);
        }
        event.metadata["pythia8:num_hadronisation_trials"] = num_trials;
      }

      if (EventInterface pythia_event(event, *pythia_); pythia_event.decay(weight)) {
        for (const auto& [role, fragments] : hadronised)  // decay/fragment rest of the event content
          pythia_event.update(role, fragments);
        return true;
      }
      return false;
    }

  private:
    void* enginePtr() override { return pythia_.get(); }

    const std::unique_ptr<Pythia8::Pythia> pythia_;  ///< Pythia 8 core to be wrapped
    const std::unique_ptr<utils::RandomGenerator> random_generator_;
    const double momentum_tolerance_;
    const double primordial_kt_;
    const std::string output_config_;

    std::unordered_map<Particle::Role, Particle::Role> roles_;  ///< (incoming, outgoing) dissociating beams
  };
}  // namespace cepgen::pythia8
using Pythia8Hadroniser = cepgen::pythia8::Hadroniser;
REGISTER_MODIFIER("pythia8", Pythia8Hadroniser);
