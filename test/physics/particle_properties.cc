/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2025  Laurent Forthomme
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

#include "CepGen/Physics/ParticleProperties.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Test.h"

using namespace std;
using namespace std::string_literals;

int main(int argc, char* argv[]) {
  cepgen::ArgumentsParser(argc, argv).parse();

  cepgen::ParametersList plist;
  plist.set<string>("name", "laurenteron"s);

  cepgen::ParticleProperties prop(plist);
  CG_TEST_EQUAL(prop.name, plist.get<string>("name"), "custom particle name");

  prop.pdgid = 42;
  CG_TEST_EQUAL(prop.parameters().get<cepgen::pdgid_t>("pdgid"), prop.pdgid, "post-defined particle id change");

  CG_TEST_SUMMARY;
}
