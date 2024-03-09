/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2024  Laurent Forthomme
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
#include "CepGen/Utils/Environment.h"
#include "CepGen/Utils/Test.h"

using namespace std;

int main(int argc, char* argv[]) {
  string path;
  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("input,i",
                           "path to the MCD file",
                           &path,
                           cepgen::utils::env::get("CEPGEN_PATH") + "/External/mass_width_2023.txt")
      .parse();

  pdg::MCDFileParser::parse(path);
  cepgen::PDG::get().dump();
  CG_TEST_SET_PRECISION(1.e-6);

  CG_TEST_EQUIV(cepgen::PDG::get().mass(cepgen::PDG::diffractiveProton), 0., "diffractive proton bare mass");
  CG_TEST_EQUIV(cepgen::PDG::get().mass(6), 172.5, "top mass");
  CG_TEST_EQUIV(cepgen::PDG::get().width(13), 2.9959836e-19, "muon width");
  CG_TEST_EQUIV(cepgen::PDG::get().mass(12), 0., "electron neutrino mass");
  CG_TEST_EQUIV(cepgen::PDG::get().mass(14), 0., "muon neutrino mass");
  CG_TEST_EQUIV(cepgen::PDG::get().mass(16), 0., "tau neutrino mass");
  {
    const auto exp_ele_ch = std::vector<double>{-1., 1.};
    CG_TEST_EQUAL(cepgen::PDG::get().charges(11), exp_ele_ch, "electron/positron charges");
  }
  CG_TEST_EQUAL(cepgen::PDG::get().charges(22), std::vector<double>{}, "photon charge");

  CG_TEST_SUMMARY;
}
