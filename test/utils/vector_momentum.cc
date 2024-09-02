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

#include "CepGen/Physics/Momentum.h"
#include "CepGen/Utils/Algebra.h"
#include "CepGen/Utils/Test.h"

using namespace std;

int main() {
  {
    const cepgen::Vector vec{1., 2., 3., 4.};
    const cepgen::Momentum mom(vec);
    CG_TEST_EQUAL(static_cast<cepgen::Vector>(mom), vec, "vec->mom->vec");
  }
  {
    const cepgen::Momentum mom{1., 2., 3., 4.};
    const cepgen::Vector vec = mom;
    CG_TEST_EQUAL(static_cast<cepgen::Momentum>(vec), mom, "mom->vec->mom");
  }
  CG_TEST_SUMMARY;
}
