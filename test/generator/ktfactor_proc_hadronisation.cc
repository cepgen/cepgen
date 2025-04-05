/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024-2025  Laurent Forthomme
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
#include "CepGen/Utils/EventUtils.h"
#include "CepGen/Utils/Test.h"

using namespace std;
using namespace std::string_literals;

int main(int argc, char* argv[]) {
  string hadroniser;
  cepgen::ArgumentsParser(argc, argv)
      .addArgument("hadroniser,H", "hadronisation/fragmentation algorithm to use", &hadroniser)
      .parse();
  auto gen = cepgen::Generator();

  auto evt = cepgen::utils::generateLPAIREvent();
  evt.dump();

  gen.runParameters().setProcess(cepgen::ProcessFactory::get().build(
      "lpair",
      cepgen::ParametersList().set(
          "kinematics",
          cepgen::ParametersList()
              .set<double>("cmEnergy", 13.e3)
              .setAs<int, cepgen::mode::Kinematics>("mode", cepgen::mode::Kinematics::InelasticElastic))));

  const auto prefix = "["s + hadroniser + "] ";

  const auto hadroniser_algo = cepgen::EventModifierFactory::get().build(hadroniser);
  CG_TEST(hadroniser_algo, prefix + "algorithm construction");
  hadroniser_algo->setCrossSection(cepgen::Value{1.46161e-1, 1.25691e-3});
  hadroniser_algo->initialise(gen.runParameters());
  double evt_weight = 1.;

  const auto evt_before_particles = evt.particles().size();
  hadroniser_algo->run(evt, evt_weight, true);
  CG_TEST(evt_weight == 1., prefix + "no event weight modification in fast mode");
  CG_TEST(evt.particles().size() == evt_before_particles, prefix + "no event modification in fast mode");

  hadroniser_algo->run(evt, evt_weight, false);

  CG_DEBUG("main") << "Hadroniser-filtered event:\n" << evt;

  CG_TEST_EQUAL(evt_weight, 1., prefix + "event weight");
  CG_TEST(evt(cepgen::Particle::Role::OutgoingBeam1).size() > 1, prefix + "decayed diffractive beam system");
  CG_TEST(evt(cepgen::Particle::Role::OutgoingBeam2).size() == 1, prefix + "undecayed elastic beam system");
  cepgen::Momentum daughters_total_momentum;
  for (const auto& daughter : evt.stableDaughters(evt(cepgen::Particle::Role::OutgoingBeam1)[0], true))
    daughters_total_momentum += daughter.get().momentum();
  CG_TEST_EQUIV((daughters_total_momentum - evt(cepgen::Particle::Role::OutgoingBeam1)[0].momentum()).p(),
                0.,
                prefix + "diffractive system momentum balance");

  CG_TEST_SUMMARY;
}
