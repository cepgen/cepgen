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

#include <cmath>
#include <functional>
#include <random>

#include "CepGen/Core/RunParameters.h"
#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventModifier.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/EventUtils.h"
#include "CepGen/Utils/Test.h"

using namespace std;
using namespace cepgen;

int main(int argc, char* argv[]) {
  CG_TEST_EXIT_GUARD();  // e.g. Fortran STOP statements exiting with a zero status
  int num_events;
  bool debug;
  ArgumentsParser(argc, argv)
      .addOptionalArgument("num-events,n", "number of events to decay", &num_events, 500)
      .addOptionalArgument("debug,d", "debugging mode", &debug, false)
      .parse();
  CG_TEST_DEBUG(debug);
  initialise();

  const RunParameters run_parameters;
  auto tauola = EventModifierFactory::get().build("tauola", ParametersList().set<int>("seed", 42));
  CG_TEST(tauola != nullptr, "Tauola interface retrieval");
  if (tauola) {
    tauola->initialise(run_parameters);

    constexpr double tolerance = 1.e-4;  // relative to the mother energy (parts of TAUOLA work in single precision)
    const auto conserves_momentum = [&](const Event& event, const Particle& mother) {
      auto sum = Momentum::fromPxPyPzE(0., 0., 0., 0.);
      for (const auto& child_id : mother.children())
        sum += event(child_id).momentum();
      const auto diff = sum - mother.momentum();
      return diff.p() + std::fabs(diff.energy()) < tolerance * mother.momentum().energy();
    };

    std::mt19937 rng(42);
    std::uniform_real_distribution<double> cos_theta(-1., 1.), phi(0., 2. * M_PI);
    size_t num_failed_runs = 0, num_undecayed = 0, num_unbalanced = 0, num_wrong_neutrino = 0, num_wrong_charge = 0,
           num_wrong_role = 0, num_modified_beams = 0;
    for (int i = 0; i < num_events; ++i) {
      auto event = utils::generateBalancedPairEvent(+(spdgid_t)PDG::tau, cos_theta(rng), phi(rng));
      const auto original_event = event;
      double weight = 1.;
      if (!tauola->run(event, weight, false)) {
        ++num_failed_runs;
        continue;
      }
      for (const auto role : {Particle::Role::OutgoingBeam1, Particle::Role::OutgoingBeam2})
        if (event.oneWithRole(role).momentum() != original_event.oneWithRole(role).momentum() ||
            !event.oneWithRole(role).children().empty())
          ++num_modified_beams;
      for (const auto& tau_ref : event[original_event.ids(Particle::Role::CentralSystem)]) {
        const auto& tau = tau_ref.get();
        if (tau.children().empty() || tau.status() != Particle::Status::Propagator) {
          ++num_undecayed;
          continue;
        }
        // tau- -> nu_tau + X-, tau+ -> anti-nu_tau + X+
        bool found_neutrino = false;
        for (const auto& child_id : tau.children())
          found_neutrino |= event(child_id).integerPdgId() == (tau.integerPdgId() > 0 ? +16 : -16);
        if (!found_neutrino)
          ++num_wrong_neutrino;
        // four-momentum and charge conservation along the whole decay chain
        double final_state_charge = 0.;
        bool charge_known = true;
        std::function<void(const Particle&)> browse = [&](const Particle& particle) {
          if (particle.role() != Particle::Role::CentralSystem)
            ++num_wrong_role;
          if (particle.children().empty()) {
            if (PDG::get().has(particle.pdgId()))
              final_state_charge += PDG::get().charge(particle.integerPdgId());
            else
              charge_known = false;
            return;
          }
          if (!conserves_momentum(event, particle))
            ++num_unbalanced;
          for (const auto& child_id : particle.children())
            browse(event(child_id));
        };
        browse(tau);
        if (charge_known && std::fabs(final_state_charge - PDG::get().charge(tau.integerPdgId())) > 1.e-6)
          ++num_wrong_charge;
      }
    }
    CG_TEST_EQUAL(num_failed_runs, 0ul, "events with a successful decay");
    CG_TEST_EQUAL(num_undecayed, 0ul, "all taus decayed");
    CG_TEST_EQUAL(num_wrong_neutrino, 0ul, "tau neutrino of the right sign in each decay");
    CG_TEST_EQUAL(num_unbalanced, 0ul, "four-momentum conservation in all decay vertices");
    CG_TEST_EQUAL(num_wrong_charge, 0ul, "charge conservation in tau decays");
    CG_TEST_EQUAL(num_wrong_role, 0ul, "decay products attached to the central system");
    CG_TEST_EQUAL(num_modified_beams, 0ul, "beam remnants untouched");

    {  // an event without tau must be left untouched, and reported as unmodified
      auto event = utils::generateBalancedPairEvent(+(spdgid_t)PDG::muon);
      const auto event_size = event.size();
      double weight = 1.;
      CG_TEST(!tauola->run(event, weight, false), "no decay reported for a tau-less event");
      CG_TEST_EQUAL(event.size(), event_size, "tau-less event untouched");
    }
  }
  CG_TEST_SUMMARY;
}
