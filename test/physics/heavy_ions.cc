/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022  Laurent Forthomme
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

#include "CepGen/Generator.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/Test.h"

int main() {
  cepgen::initialise();
  const auto mp = cepgen::PDG::get().mass(2212), mn = cepgen::PDG::get().mass(2112);
  {
    const cepgen::HeavyIon proton(2212);
    CG_TEST_EQUAL(proton.mass(), mp, "single proton mass");
    CG_TEST_EQUAL((cepgen::pdgid_t)proton, 2212, "single proton PDG id");
  }
  {
    const cepgen::HeavyIon neutron(2112);
    CG_TEST_EQUAL(neutron.mass(), mn, "single neutron mass");
    CG_TEST_EQUAL((cepgen::pdgid_t)neutron, 2112, "single neutron PDG id");
  }
  {
    const auto hi = cepgen::HeavyIon::Pb();
    CG_TEST_EQUAL(hi.massP(), (int)hi.Z * mp, "proton masses in lead ion");
    CG_TEST_EQUAL(hi.massN(), (hi.A - (int)hi.Z) * mn, "neutron masses in lead ion");
  }
  CG_TEST_SUMMARY;
}
