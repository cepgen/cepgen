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

#include "CepGen/Core/Exception.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Test.h"
#include "CepGenMadGraph/Process.h"

using namespace std;
using namespace std::string_literals;

int main(int argc, char* argv[]) {
  cepgen::ArgumentsParser(argc, argv).parse();
  cepgen::initialise();

  const cepgen::pdgid_t my_part = 13;
  const auto my_part_mass = 42.;

  cepgen::PDG::get().define(
      cepgen::ParticleProperties(my_part, "la", "laurentino", 0., my_part_mass, 0., {-3, 3}, true));

  auto mg5 = cepgen::ProcessFactory::get().build(
      "mg5_aMC",
      cepgen::ParametersList()
          .set("kinematicsGenerator", cepgen::ParametersList().setName("coll:2to4"s))
          .set("extraParticles", cepgen::ParametersList().set("la", cepgen::PDG::get()(my_part)))
          .set("process", "a a > la+ la-"s));
  mg5->initialise();

  const auto& proc_evt = mg5->event();
  CG_TEST_EQUAL(proc_evt.oneWithRole(cepgen::Particle::Role::Parton1).pdgId(), cepgen::PDG::photon, "parton 1 PDG id");
  CG_TEST_EQUAL(proc_evt.oneWithRole(cepgen::Particle::Role::Parton2).pdgId(), cepgen::PDG::photon, "parton 2 PDG id");
  CG_TEST_EQUAL(proc_evt(cepgen::Particle::Role::CentralSystem).size(), 2, "cent.part.multiplicity");
  for (size_t i = 0; i < 2; ++i) {
    const auto& cent = proc_evt(cepgen::Particle::Role::CentralSystem).at(i);
    CG_TEST_EQUAL(cent.integerPdgId(), short((i == 0 ? -1 : +1) * my_part), "cent." + to_string(i) + " PDG id");
    CG_TEST_EQUAL(cent.momentum().mass(), my_part_mass, "cent." + to_string(i) + " mass");
  }

  CG_TEST_SUMMARY;
}
