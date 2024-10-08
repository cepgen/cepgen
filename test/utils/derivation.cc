/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2021-2022  Laurent Forthomme
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

#include "CepGen/Generator.h"
#include "CepGen/Modules/DerivatorFactory.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Derivator.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/Test.h"

using namespace std;

int main(int argc, char* argv[]) {
  string plotter;
  vector<string> derivators;

  cepgen::initialise();

  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("plotter,p", "type of plotter to use", &plotter, "text")
      .addOptionalArgument(
          "derivators,D", "type of derivators to use", &derivators, cepgen::DerivatorFactory::get().modules())
      .parse();

  auto plt = cepgen::DrawerFactory::get().build(plotter);
  for (const auto& derivator_name : derivators) {
    auto der = cepgen::DerivatorFactory::get().build(derivator_name, cepgen::ParametersList().set<double>("h", 0.05));

    // test 1D graph
    cepgen::utils::Graph1D graph_sin("graph_sin", "sin(x)"), graph_cos("graph_cos", "cos(x)"),
        graph_der_sin("graph_der_sin", "(sin(x))'"), graph_diff("graph_diff", "cos(x)-(sin(x))'");
    for (const auto x : cepgen::Limits{-M_PI, M_PI}.generate(25)) {
      graph_sin.addPoint(x, sin(x));
      graph_cos.addPoint(x, cos(x));
      const auto der_sin = der->derivate([](double x) { return sin(x); }, x);
      graph_der_sin.addPoint(x, der_sin);
      graph_diff.addPoint(x, cos(x) - der_sin);
    }
    (void)plt->draw({&graph_sin, &graph_der_sin, &graph_diff}, "test_deriv_" + derivator_name);

    const auto chi2 = graph_cos.chi2(graph_der_sin);
    CG_TEST(chi2 <= 1.e-6, "chi^2 test for " + derivator_name);
  }

  CG_TEST_SUMMARY;
}
