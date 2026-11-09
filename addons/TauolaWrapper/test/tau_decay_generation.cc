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

#include "CepGen/Core/RunParameters.h"
#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventModifier.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Modules/PhaseSpaceGeneratorFactory.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Collections.h"
#include "CepGen/Utils/Test.h"

using namespace std;
using namespace std::string_literals;
using namespace cepgen;

/// Generate tau pairs with CepGen processes, decay them with Tauola++, and check the resulting events
int main(int argc, char* argv[]) {
  CG_TEST_EXIT_GUARD();  // e.g. Fortran STOP statements exiting with a zero status
  vector<string> processes;
  int num_events, num_calls;
  double pt_min, eta_max;
  bool debug;
  ArgumentsParser(argc, argv)
      .addOptionalArgument("processes,p", "processes to generate", &processes, vector<string>{"lpair", "pptoff"})
      .addOptionalArgument("num-events,n", "number of events to generate per process", &num_events, 100)
      .addOptionalArgument("num-calls,c", "number of function calls per integration iteration", &num_calls, 50000)
      .addOptionalArgument("ptmin", "minimum transverse momentum of the central particles", &pt_min, 5.)
      .addOptionalArgument("etamax", "maximum absolute pseudorapidity of the central particles", &eta_max, 2.5)
      .addOptionalArgument("debug,d", "debugging mode", &debug, false)
      .parse();
  CG_TEST_DEBUG(debug);
  initialise();

  // generator, steered for a fast integration of a 13.6 TeV elastic gamma-gamma -> X+X- process
  Generator gen;
  if (auto& integrator = gen.runParameters().integrator(); !integrator.hasName() || integrator.name() == "Vegas")
    integrator.setName("Vegas").set<int>("numFunctionCalls", num_calls);
  const auto configure_run = [&](const string& process_name) {
    auto process_parameters = ParametersList().set<int>("pair", (int)PDG::tau);
    if (process_name != "lpair")  // kt-factorised processes
      process_parameters.set("kinematicsGenerator",
                             PhaseSpaceGeneratorFactory::get().describeParameters("kt:2to4").parameters());
    gen.runParameters().setProcess(ProcessFactory::get().build(process_name, process_parameters));
    gen.runParameters().process().kinematics().setParameters(ParametersList()
                                                                 .set<vector<int> >("pdgIds", {2212, 2212})
                                                                 .set<double>("sqrtS", 13.6e3)
                                                                 .set<int>("mode", 1 /* elastic-elastic */)
                                                                 .set<double>("ptmin", pt_min)
                                                                 .set<Limits>("eta", {-eta_max, eta_max}));
  };

  // reference cross sections, computed without any event modification algorithm
  unordered_map<string, Value> reference_cross_sections;
  for (const auto& process_name : processes) {
    configure_run(process_name);
    reference_cross_sections[process_name] = gen.computeXsection();
  }

  auto tauola = EventModifierFactory::get().build("tauola", ParametersList().set<int>("seed", 42));
  CG_TEST(tauola != nullptr, "Tauola interface retrieval");
  if (tauola) {
    gen.runParameters().addModifier(std::move(tauola));

    for (const auto& process_name : processes) {
      const auto prefix = "["s + process_name + "] ";
      configure_run(process_name);
      const auto cross_section = gen.computeXsection();
      // tau decays have a unit weight, and the cuts are only meant to act on the taus themselves
      CG_TEST_VALUES(cross_section, reference_cross_sections.at(process_name), 3., prefix + "cross section unchanged");
      if (cross_section < 0.5 * reference_cross_sections.at(process_name)) {
        CG_WARNING("main") << prefix << "Cross section collapsed after tau decays; skipping the events generation.";
        continue;
      }

      size_t num_generated = 0, num_taus = 0, num_undecayed = 0, num_wrong_neutrino = 0, num_unbalanced = 0,
             num_products_outside_cuts = 0;
      gen.generate(num_events, [&](const Event& event, size_t) {
        ++num_generated;
        const auto intermediate_id = event(Particle::Role::Intermediate).at(0).id();
        for (const auto& tau : event(Particle::Role::CentralSystem)) {
          if (tau.pdgId() != PDG::tau || tau.mothers().count(intermediate_id) == 0)
            continue;  // only consider the taus produced by the process
          ++num_taus;
          if (tau.children().empty()) {
            ++num_undecayed;
            continue;
          }
          bool found_neutrino = false;
          for (const auto& child_id : tau.children())
            found_neutrino |= event(child_id).integerPdgId() == (tau.integerPdgId() > 0 ? +16 : -16);
          if (!found_neutrino)
            ++num_wrong_neutrino;
          function<void(const Particle&)> browse = [&](const Particle& particle) {
            if (particle.children().empty()) {  // stable decay product
              if (particle.momentum().pt() < pt_min || std::fabs(particle.momentum().eta()) > eta_max)
                ++num_products_outside_cuts;
              return;
            }
            Momentum sum;
            for (const auto& child_id : particle.children())
              sum += event(child_id).momentum();
            if ((sum - particle.momentum()).p() + std::fabs(sum.energy() - particle.momentum().energy()) >
                1.e-4 * particle.momentum().energy())
              ++num_unbalanced;
            for (const auto& child_id : particle.children())
              browse(event(child_id));
          };
          browse(tau);
        }
      });
      CG_TEST_EQUAL(num_generated, (size_t)num_events, prefix + "number of events generated");
      CG_TEST_EQUAL(num_taus, 2 * num_generated, prefix + "two taus produced per event");
      CG_TEST_EQUAL(num_undecayed, 0ul, prefix + "all taus decayed");
      CG_TEST_EQUAL(num_wrong_neutrino, 0ul, prefix + "tau neutrino of the right sign in each decay");
      CG_TEST_EQUAL(num_unbalanced, 0ul, prefix + "four-momentum conservation in all decay vertices");
      CG_TEST(num_products_outside_cuts > 0, prefix + "decay products not subject to the central system cuts");
    }
  }
  CG_TEST_SUMMARY;
}
