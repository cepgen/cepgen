/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2022  Laurent Forthomme
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

#include "CepGen/Physics/MCDFileParser.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Test.h"

using namespace std;

int main(int argc, char* argv[]) {
  string path;
  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("input,i", "path to the MCD file", &path, "../External/mass_width_2023.txt")
      .parse();

  pdg::MCDFileParser::parse(path);
  cepgen::PDG::get().dump();

  CG_TEST_EQUAL(cepgen::PDG::get().mass(cepgen::PDG::diffractiveProton), 0., "diffractive proton bare mass");
  CG_TEST_EQUAL(cepgen::PDG::get().mass(6), 172.5, "top mass");
  CG_TEST_EQUAL(cepgen::PDG::get().width(13), 2.9959836e-19, "muon width");

  CG_TEST_SUMMARY;
}
