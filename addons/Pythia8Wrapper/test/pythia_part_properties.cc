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
#include "CepGenPythia8/EventInterface.h"

using namespace std;

int main(int argc, char* argv[]) {
  vector<int> pdgids;
  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("particles,p", "list of PDGids to probe", &pdgids, vector<int>{})
      .parse();
  //cepgen::initialise();  // do not initialise to fetch attributes from Pythia directly!
  auto pythia8 = make_unique<Pythia8::Pythia>();
  if (pdgids.empty())
    for (const auto& pd : pythia8->particleData)
      if (pd.first != 0)
        pdgids.emplace_back(pd.first);
  for (const auto& pdgid : pdgids) {
    {
      if (cepgen::PDG::get().has(pdgid))  // skip particles already defined natively by CepGen
        continue;
      auto unknown_pdgid = [&pdgid]() { cepgen::PDG::get()(pdgid); };
      CG_TEST_EXCEPT(unknown_pdgid, "unknown PDG id [" + to_string(pdgid) + "]");
    }

    if (auto data = pythia8->particleData.findParticle(pdgid); data) {
      cepgen::pythia8::EventInterface::checkPDGid(*data);
      const auto name = cepgen::PDG::get().name(pdgid);
      CG_TEST(!name.empty(), "valid name [" + to_string(pdgid) + "=" + name + "]");
      const auto mass = cepgen::PDG::get().mass(pdgid);
      CG_TEST(mass >= 0., "valid mass [" + to_string(pdgid) + "=" + to_string(mass) + "]");
    }
  }

  CG_TEST_SUMMARY;
}
