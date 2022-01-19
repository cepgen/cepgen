/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2022  Laurent Forthomme
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
#include <iostream>
#include <random>

#include "CepGen/Core/Exception.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/GSLDerivator.h"
#include "CepGen/Utils/Graph.h"

using namespace std;

int main(int argc, char* argv[]) {
  string plotter;

  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("plotter,p", "type of plotter to user", &plotter, "text")
      .parse();
  cepgen::initialise();

  auto plt = cepgen::utils::DrawerFactory::get().build(plotter);
  auto der = cepgen::utils::GSLDerivator(cepgen::ParametersList().set<double>("h", 0.05));

  // test 1D graph
  cepgen::utils::Graph1D graph_sin("graph_sin"), graph_der_sin("graph_der_sin"), graph_diff("graph_diff");
  for (double x = -M_PI; x <= M_PI; x += 0.25) {
    graph_sin.addPoint(x, sin(x));
    const auto der_sin = der.eval([](double x) { return sin(x); }, x);
    graph_der_sin.addPoint(x, der_sin);
    graph_diff.addPoint(x, cos(x) - der_sin);
  }
  plt->draw({&graph_sin, &graph_der_sin, &graph_diff}, "test_deriv");

  return 0;
}
