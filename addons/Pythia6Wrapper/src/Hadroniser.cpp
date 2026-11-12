/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2026  Laurent Forthomme
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
#include <CepGen/Utils/String.h>

#include "CepGenPythia6/DiffractiveState.h"
#include "CepGenPythia6/EventInterface.h"
#include "CepGenPythia6/Pythia6Interface.h"

namespace cepgen::pythia6 {
  /// Interface to the Pythia 6 algorithm
  /// \note It can be used in a single particle decay mode as well as a full event hadronisation using the string model, as in Jetset.
  class Hadroniser final : public hadr::Hadroniser {
  public:
    explicit Hadroniser(const ParametersList& params)
        : hadr::Hadroniser(params),
          random_generator_(RandomGeneratorFactory::get().build(steer<ParametersList>("randomGenerator"))),
          momentum_tolerance_(steer<double>("momentumTolerance")),
          primordial_kt_(steer<double>("primordialKt")) {}

    static ParametersDescription description() {
      auto desc = hadr::Hadroniser::description();
      desc.setDescription("Interface to the Pythia 6 string hadronisation/fragmentation algorithm");
      desc.add("maxTrials", 10).setDescription("maximum number of hadronisation attempts per dissociated system");
      desc.add("randomGenerator", RandomGeneratorFactory::get().describeParameters("stl"))
          .setDescription("random number generator to use for the partonic content of the dissociated systems");
      desc.add("momentumTolerance", 1.e-6)
          .setDescription("maximum fractional momentum tolerance to consider a hadronisation successful");
      desc.add("primordialKt", 0.)
          .setDescription("RMS of the spectator transverse momentum relative to the beam axis (GeV)");
      return desc;
    }

    void readString(const std::string& param) override { pygive(param); }
    void initialise() override {
      CG_WARNING("pythia6:Hadroniser") << "Branching fraction not yet implemented in this hadroniser.\n\t"
                                       << "You will have to specify manually the multiplication factor according\n\t"
                                       << "to your list of open channels.";
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
      std::vector<std::pair<Particle::Role, hadr::Fragments> > hadronised;
      if (!fast && fragment_remnants_) {  // do not treat beam remnants when running in fast mode
        int num_trials = 0;
        for (const auto& [incoming, outgoing] : roles_) {
          if (!event.hasRole(outgoing) || event(outgoing).size() != 1)  // already fragmented by another algorithm
            continue;
          const pythia6::DiffractiveState diffractive_state(
              *random_generator_, event.oneWithRole(incoming), event.oneWithRole(outgoing).momentum(), primordial_kt_);
          if (!diffractive_state.valid()) {
            CG_WARNING("pythia6:Hadroniser") << "Invalid remnants kinematics for " << outgoing << " (mass "
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
          hadronised.emplace_back(outgoing, std::move(fragments));
        }
        event.metadata["pythia6:num_hadronisation_trials"] = num_trials;
      }

      // then decay/fragment the rest of the event content (central system, strings)
      EventInterface pythia_event(event);
      if (!pythia_event.decay())
        return false;
      CG_DEBUG_LOOP("pythia6:Hadroniser") << utils::s("string object", pythia_event.numStrings(), true)
                                          << " identified and constructed in the event content.";

      for (const auto& [role, fragments] : hadronised)
        pythia_event.update(role, fragments);
      return true;
    }

  private:
    const std::unique_ptr<utils::RandomGenerator> random_generator_;
    const double momentum_tolerance_;
    const double primordial_kt_;

    std::unordered_map<Particle::Role, Particle::Role> roles_;  ///< (incoming, outgoing) dissociating beams
  };
}  // namespace cepgen::pythia6
using Pythia6Hadroniser = cepgen::pythia6::Hadroniser;
REGISTER_MODIFIER("pythia6", Pythia6Hadroniser);
