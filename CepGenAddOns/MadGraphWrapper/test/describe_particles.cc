/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024  Laurent Forthomme
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
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Test.h"
#include "CepGenAddOns/MadGraphWrapper/Utils.h"

using namespace std;
using namespace std::string_literals;

int main(int argc, char* argv[]) {
  cepgen::ArgumentsParser(argc, argv).parse();
  cepgen::initialise();

  const cepgen::pdgid_t my_part = 13;
  {
    const auto pprop = cepgen::mg5amc::describeParticle("a", "sm");
    CG_LOG << "photon:" << pprop;
  }
  for (const auto& part : vector<pair<cepgen::pdgid_t, string> >{{11, "e"s}, {13, "mu"s}, {15, "ta"s}}) {
    const auto pprop = cepgen::mg5amc::describeParticle(part.second, "sm");
    const auto cprop = cepgen::PDG::get()(part.first);
    CG_TEST_EQUAL(pprop.pdgid, part.first, part.second + " PDG");
    CG_TEST_EQUAL(pprop.fermion, cprop.fermion, part.second + " fermion");
    CG_TEST_EQUAL(pprop.mass, cprop.mass, part.second + " mass");
  }
  //CG_TEST_EQUAL(cent.momentum().mass(), my_part_mass, "cent." + to_string(i) + " mass");

  CG_TEST_SUMMARY;
}
