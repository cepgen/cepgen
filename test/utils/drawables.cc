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

#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/Histogram.h"
#include "CepGen/Utils/Test.h"

using namespace std;

int main() {
  auto valid_hist1d = cepgen::utils::Hist1D(10, {0., 1.});
  CG_TEST_EQUAL(valid_hist1d.nbins(), 10, "number of bins in 1D histogram");

  valid_hist1d.fill(2.);
  CG_TEST_EQUAL(valid_hist1d.overflow(), 1, "overflow counter");
  CG_TEST_EQUAL(valid_hist1d.underflow(), 0, "underflow counter");
  CG_TEST_EQUAL(valid_hist1d.integral(), 0, "integral");
  CG_TEST_EQUAL(valid_hist1d.integral(true), 1, "integral (with out-of-range)");

  valid_hist1d.add(valid_hist1d, 2);
  CG_TEST_EQUAL(valid_hist1d.integral(), 0, "integral (tripled 1D histogram)");
  CG_TEST_EQUAL(valid_hist1d.integral(true), 3, "integral (with out-of-range, tripled 1D histogram)");

  auto valid_hist2d = cepgen::utils::Hist2D(10, {0., 1.}, 20, {0., 1.});
  CG_TEST_EQUAL(valid_hist2d.nbinsX(), 10, "number of x-bins in 2D histogram");
  CG_TEST_EQUAL(valid_hist2d.nbinsY(), 20, "number of y-bins in 2D histogram");

  valid_hist2d.fill(2., 2.);
  valid_hist2d.fill(-2., -2.);
  valid_hist2d.add(valid_hist2d, 2);
  CG_TEST_EQUAL(valid_hist2d.integral(), 0, "integral (tripled 2D histogram)");
  CG_TEST_EQUAL(valid_hist2d.integral(true), 6, "integral (with out-of-range, tripled 2D histogram)");

  {
    auto invalid_hist1d = []() { cepgen::utils::Hist1D empty_hist(0, {0., 1.}); };
    CG_TEST_EXCEPT(invalid_hist1d, "zero-binned 1D histogram");
    auto invalid_hist2d = []() { cepgen::utils::Hist2D empty_hist(1, {0., 1.}, 0, {0., 1.}); };
    CG_TEST_EXCEPT(invalid_hist2d, "zero-binned 2D histogram");
  }

  cepgen::utils::Graph1D empty_graph("empty graph");

  CG_TEST_SUMMARY;
}
