/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2024  Laurent Forthomme
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

#include "CepGen/Core/RunParameters.h"
#include "CepGen/EventFilter/EventModifier.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Physics/Beam.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Test.h"
#include "CepGenAddOns/Common/EventUtils.h"

using namespace std;

int main(int argc, char* argv[]) {
  int seed;
  cepgen::ArgumentsParser(argc, argv).addOptionalArgument("seed,s", "hadroniser seed", &seed, -1).parse();
  auto gen = cepgen::Generator();

  auto evt = cepgen::utils::generateLPAIREvent();
  evt.dump();

  gen.runParameters().setProcess(cepgen::ProcessFactory::get().build(
      "lpair",
      cepgen::ParametersList().set(
          "kinematics",
          cepgen::ParametersList()
              .set<double>("cmEnergy", 13.e3)
              .set<vector<int> >("pdgIds", {2212, 2212})
              .setAs<int, cepgen::mode::Kinematics>("mode", cepgen::mode::Kinematics::InelasticElastic))));

  if (seed < 0)
    seed = std::time(0) % 100'000;

  auto cg_pythia =
      cepgen::EventModifierFactory::get().build("pythia8", cepgen::ParametersList().set<int>("seed", seed));
  cg_pythia->setCrossSection(cepgen::Value{1.46161e-1, 1.25691e-3});
  cg_pythia->initialise(gen.runParameters());
  double evt_weight = 1.;

  const auto evt_before_particles = evt.particles().size();
  cg_pythia->run(evt, evt_weight, true);
  CG_TEST_EQUAL(evt_weight, 1., "no event weight modification in fast mode");
  CG_TEST(evt.particles().size() == evt_before_particles, "no event modification in fast mode");

  cg_pythia->run(evt, evt_weight, false);

  CG_DEBUG("main") << "Pythia 8-filtered event:\n" << evt;

  CG_TEST_EQUAL(evt_weight, 1., "event weight");
  CG_TEST(evt.stableParticlesWithRole(cepgen::Particle::Role::OutgoingBeam1).size() > 1,
          "decayed diffractive beam system");
  CG_TEST(evt.stableParticlesWithRole(cepgen::Particle::Role::OutgoingBeam2).size() == 1,
          "undecayed elastic beam system");
  cepgen::Momentum daugh_total_momentum;
  for (const auto& daugh : evt.stableDaughters(evt(cepgen::Particle::Role::OutgoingBeam1)[0], true))
    daugh_total_momentum += daugh.get().momentum();
  CG_LOG << daugh_total_momentum << ":" << evt(cepgen::Particle::Role::OutgoingBeam1)[0].momentum();
  CG_TEST_EQUIV((daugh_total_momentum - evt(cepgen::Particle::Role::OutgoingBeam1)[0].momentum()).p(),
                0.,
                "diffractive system momentum balance");

  CG_TEST_SUMMARY;
}
