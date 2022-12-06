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

#include "CepGen/Core/Exception.h"
#include "CepGen/Generator.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/Message.h"

int main() {
  cepgen::initialise();
  const auto mp = cepgen::PDG::get().mass(2212), mn = cepgen::PDG::get().mass(2112);
  {
    const cepgen::HeavyIon proton(2212);
    if (proton.mass() != mp)
      throw CG_FATAL("main") << "Test 1.1 failed!";
    if ((cepgen::pdgid_t)proton != 2212)
      throw CG_FATAL("main") << "Test 1.2 failed!";
  }
  {
    const cepgen::HeavyIon neutron(2112);
    if (neutron.mass() != mn)
      throw CG_FATAL("main") << "Test 2.1 failed!";
    if ((cepgen::pdgid_t)neutron != 2112)
      throw CG_FATAL("main") << "Test 2.2 failed!";
  }
  {
    const auto hi = cepgen::HeavyIon::Pb();
    if (hi.massP() != (int)hi.Z * mp)
      throw CG_FATAL("main") << "Test 3.1 failed!";
    if (hi.massN() != (hi.A - (int)hi.Z) * mn)
      throw CG_FATAL("main") << "Test 3.2 failed!";
  }
  return 0;
}
