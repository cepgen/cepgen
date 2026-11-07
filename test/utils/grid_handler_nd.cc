/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2026  Laurent Forthomme
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

#include <cmath>
#include <random>
#include <string>
#include <vector>

#include "CepGen/Utils/GridHandler.h"
#include "CepGen/Utils/Test.h"

using namespace std;

/// A multilinear function (linear in each coordinate separately), exactly reproduced by a multilinear interpolation
template <size_t D>
double multilinear(const vector<double>& x) {
  double out = 0.5, prod = 1.;
  for (size_t i = 0; i < D; ++i) {
    out += (0.1 + 0.05 * i) * x[i];
    prod *= x[i];
    if (i > 0)
      out -= 0.02 * x[i - 1] * x[i];
  }
  return out + 0.01 * prod;
}

/// Fill a D-dimensional grid with non-uniformly spaced nodes, and check its interpolation inside and outside its range
template <size_t D>
void testGrid(mt19937& rng) {
  const vector<double> nodes{0., 0.1, 0.35, 0.5, 0.9, 1.};  // non-uniform spacing
  cepgen::GridHandler<D, 1> grid(cepgen::GridType::linear);
  vector<size_t> idx(D, 0);
  while (true) {
    vector<double> coord(D);
    for (size_t i = 0; i < D; ++i)
      coord[i] = nodes[idx[i]];
    grid.insert(coord, {multilinear<D>(coord)});
    size_t d = 0;
    while (d < D && ++idx[d] == nodes.size())
      idx[d++] = 0;
    if (d == D)
      break;
  }
  grid.initialise();
  const auto dim = to_string(D) + "D grid";
  CG_TEST_EQUAL(grid.values().size(), static_cast<size_t>(pow(nodes.size(), D)), dim + ": number of nodes");
  uniform_real_distribution<double> inside(0., 1.), outside(-0.5, 1.5);
  for (size_t it = 0; it < 20; ++it) {
    vector<double> coord(D);
    for (auto& c : coord)
      c = inside(rng);
    CG_TEST_EQUIV(grid.eval(coord).at(0), multilinear<D>(coord), dim + ": interpolation within range");
  }
  for (size_t it = 0; it < 10; ++it) {
    vector<double> coord(D), clamped(D);
    for (size_t i = 0; i < D; ++i) {
      coord[i] = outside(rng);
      clamped[i] = min(1., max(0., coord[i]));
    }
    CG_TEST_EQUIV(grid.eval(coord).at(0), multilinear<D>(clamped), dim + ": closest boundary value out of range");
  }
}

int main() {
  mt19937 rng(42);
  CG_TEST_SET_PRECISION(1.e-10);
  testGrid<3>(rng);
  testGrid<4>(rng);
  testGrid<5>(rng);
  CG_TEST_SUMMARY;
}
