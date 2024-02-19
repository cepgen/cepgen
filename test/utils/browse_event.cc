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

#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventBrowser.h"
#include "CepGen/Generator.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Test.h"
#include "CepGenAddOns/Common/EventUtils.h"

using namespace std;

int main(int argc, char* argv[]) {
  cepgen::initialise();
  cepgen::ArgumentsParser(argc, argv).parse();

  cepgen::utils::EventBrowser bws;

  const auto evt = cepgen::utils::generateLPAIREvent();
  CG_LOG << evt;

  vector<pair<string, double> > values = {
      {"pdg(ib1)", evt.oneWithRole(cepgen::Particle::Role::IncomingBeam1).integerPdgId()},
      {"m(4)", evt(4).momentum().mass()},
      {"m2(4)", evt(4).momentum().mass2()},
      {"m(ob1)", evt.oneWithRole(cepgen::Particle::Role::OutgoingBeam1).momentum().mass()},
      {"acop(7,8)", 1. - fabs(evt(7).momentum().deltaPhi(evt(8).momentum())) * M_1_PI},
      {"m(7,8)", evt(4).momentum().mass()}};
  for (const auto& val_pair : values)
    CG_TEST_EQUIV(bws.get(evt, val_pair.first), val_pair.second, val_pair.first);
  CG_TEST_SUMMARY;
}
