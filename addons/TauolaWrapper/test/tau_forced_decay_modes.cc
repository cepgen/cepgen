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

#include <algorithm>
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
      .addOptionalArgument("num-events,n", "number of events to decay", &num_events, 200)
      .addOptionalArgument("debug,d", "debugging mode", &debug, false)
      .parse();
  CG_TEST_DEBUG(debug);
  initialise();

  // force tau- -> e- anti-nu_e nu_tau and tau+ -> mu+ nu_mu anti-nu_tau (no radiation by default)
  const RunParameters run_parameters;
  auto tauola = EventModifierFactory::get().build("tauola",
                                                  ParametersList()
                                                      .set<int>("seed", 4242)
                                                      .set<int>("decayingParticle", +(int)PDG::tau)
                                                      .set<int>("sameParticleDecayMode", 1 /* electron */)
                                                      .set<int>("oppositeParticleDecayMode", 2 /* muon */));
  CG_TEST(tauola != nullptr, "Tauola interface retrieval");
  if (tauola) {
    tauola->initialise(run_parameters);

    const auto decay_products = [](const Event& event, const Particle& particle) {
      std::vector<long> pdg_ids;
      for (const auto& child_id : particle.children())
        pdg_ids.emplace_back(event(child_id).integerPdgId());
      std::sort(pdg_ids.begin(), pdg_ids.end());
      return pdg_ids;
    };
    const std::vector<long> expected_tau_minus{-12, +11, +16}, expected_tau_plus{-16, -13, +14};

    std::mt19937 rng(4242);
    std::uniform_real_distribution<double> cos_theta(-1., 1.), phi(0., 2. * M_PI);
    size_t num_failed_runs = 0, num_wrong_tau_minus = 0, num_wrong_tau_plus = 0;
    for (int i = 0; i < num_events; ++i) {
      auto event = utils::generateBalancedPairEvent(+static_cast<spdgid_t>(PDG::tau), cos_theta(rng), phi(rng));
      const auto tau_ids = event.ids(Particle::Role::CentralSystem);
      double weight = 1.;
      if (!tauola->run(event, weight, false)) {
        ++num_failed_runs;
        continue;
      }
      for (const auto& tau_ref : event[tau_ids]) {
        const auto& tau = tau_ref.get();
        const auto products = decay_products(event, tau);
        if (tau.integerPdgId() > 0 && products != expected_tau_minus)
          ++num_wrong_tau_minus;
        else if (tau.integerPdgId() < 0 && products != expected_tau_plus)
          ++num_wrong_tau_plus;
      }
    }
    CG_TEST_EQUAL(num_failed_runs, 0ul, "events with a successful decay");
    CG_TEST_EQUAL(num_wrong_tau_minus, 0ul, "tau- forced into the electron mode");
    CG_TEST_EQUAL(num_wrong_tau_plus, 0ul, "tau+ forced into the muon mode");
  }
  CG_TEST_SUMMARY;
}
