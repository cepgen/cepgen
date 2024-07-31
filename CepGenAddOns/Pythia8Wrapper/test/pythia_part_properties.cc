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

#include <Pythia8/Pythia.h>

#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Test.h"
#include "CepGenAddOns/Pythia8Wrapper/EventInterface.h"

using namespace std;

int main(int argc, char* argv[]) {
  vector<int> pdgids;
  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("particles,p", "list of PDGids to probe", &pdgids, vector<int>{321, 13, 17, 42})
      .parse();
  //cepgen::initialise();  // do not initialise to fetch attributes from Pythia directly!
  auto pythia8 = make_unique<Pythia8::Pythia>();
  for (const auto& pdgid : pdgids) {
    {
      auto unknown_pdgid = [&pdgid]() { cepgen::PDG::get()(pdgid); };
      CG_TEST_EXCEPT(unknown_pdgid, "unknown PDG id [" + to_string(pdgid) + "]");
    }

    if (auto data = pythia8->particleData.findParticle(pdgid); data) {
      cepgen::pythia8::EventInterface::checkPDGid(*data);
      CG_TEST(!cepgen::PDG::get().name(pdgid).empty(), "valid name [" + to_string(pdgid) + "]");
      CG_TEST(cepgen::PDG::get().mass(pdgid) >= 0., "valid mass [" + to_string(pdgid) + "]");
    }
  }

  CG_TEST_SUMMARY;
}
