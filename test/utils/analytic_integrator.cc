/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2023  Laurent Forthomme
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
#include "CepGen/Integration/AnalyticIntegrator.h"
#include "CepGen/Modules/AnalyticIntegratorFactory.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/Test.h"

using namespace std;

int main(int argc, char* argv[]) {
  string integrator, plotter;
  int npoints;

  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("integrator,i", "analytical integrator to use", &integrator, "gsl")
      .addOptionalArgument("plotter,p", "type of plotter to user", &plotter, "")
      .addOptionalArgument("npoints,n", "number of points to compute", &npoints, 100)
      .parse();
  cepgen::initialise();

  auto integ = cepgen::AnalyticIntegratorFactory::get().build(integrator);

  // test 1D graph
  cepgen::utils::Graph1D graph_sin("graph_sin", "sin(x)"), graph_cos("graph_cos", "cos(x)"),
      graph_int_cos("graph_int_cos", "\\int_{0}^{\\pi}(cos(x))"),
      graph_diff("graph_diff", "sin(x)-\\int_{0}^{\\pi}(cos(x))'");
  for (const auto& x : cepgen::Limits{1.e-4, 2. * M_PI}.generate(npoints)) {
    graph_sin.addPoint(x, sin(x));
    graph_cos.addPoint(x, cos(x));
    const auto int_cos = integ->integrate([](double x) { return cos(x); }, cepgen::Limits{0., x});
    graph_int_cos.addPoint(x, int_cos);
    graph_diff.addPoint(x, sin(x) - int_cos);
  }
  if (!plotter.empty()) {
    auto plt = cepgen::DrawerFactory::get().build(plotter);
    (void)plt->draw({&graph_sin, &graph_int_cos, &graph_diff}, "test_deriv");
  }

  const auto chi2 = graph_sin.chi2(graph_int_cos);
  CG_TEST(chi2 <= 1.e-6, "chi^2 test");

  CG_TEST_SUMMARY;
}
