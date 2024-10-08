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
#include <random>

#include "CepGen/Generator.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/Histogram.h"
#include "CepGen/Utils/Test.h"

using namespace std;

int main(int argc, char* argv[]) {
  vector<string> plotters;

  cepgen::initialise();
  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("plotters,p", "type of plotter to user", &plotters, cepgen::DrawerFactory::get().modules())
      .parse();

  for (const auto& plotter : plotters) {
    auto plt = cepgen::DrawerFactory::get().build(plotter);

    CG_LOG << "---------- 1D graph ----------";

    // test 1D graph
    cepgen::utils::Graph1D graph1d("graph1d", "A graph of the sin(x) function");
    for (const auto x : cepgen::Limits{-M_PI, M_PI}.generate(25))
      graph1d.addPoint(x, sin(x));
    graph1d.xAxis().setLabel("x (rad)");
    graph1d.yAxis().setLabel("sin(x)");
    (void)plt->draw(graph1d);

    CG_LOG << "\n---------- 2D graph ----------";

    // test 2D graph
    cepgen::utils::Graph2D graph2d("graph2d");
    for (const auto x : cepgen::Limits{-5., 5.}.generate(21))
      for (const auto y : cepgen::Limits{-5., 5.}.generate(51))
        graph2d.addPoint(x, y, (sin(x) / x) * (sin(y) / y));
    (void)plt->draw(graph2d);

    default_random_engine gen;

    CG_LOG << "\n-------- 1D histogram --------";

    // test 1D histogram
    cepgen::utils::Hist1D hist1d(20, {-5., 5.}, "hist1d");
    cauchy_distribution<double> bw(0., 1.);
    for (size_t i = 0; i < 10000; ++i)
      hist1d.fill(bw(gen));
    hist1d.xAxis().setLabel("Random variable");
    hist1d.yAxis().setLabel("Occurrences");
    (void)plt->draw(hist1d, cepgen::utils::Drawer::Mode::logy);

    CG_LOG << "\n-------- 2D histogram --------";

    // test 2d histogram
    cepgen::utils::Hist2D hist2d(20, {-5., 5.}, 50, {-5., 5.}, "hist2d", "$\\sqrt{s} = 14$ TeV");
    normal_distribution<double> gaussian1(0., 1.), gaussian2(0., 1.);
    for (size_t i = 0; i < 1000; ++i)
      for (size_t j = 0; j < 1000; ++j)
        hist2d.fill(gaussian1(gen), gaussian2(gen));
    hist2d.xAxis().setLabel("$4\\pi\\alpha_{EM}$");
    hist2d.yAxis().setLabel("$\\Sigma(1\\pm\\epsilon)$");
    (void)plt->draw(hist2d, cepgen::utils::Drawer::Mode::logz);

    CG_LOG << "\n--------- multiplots ---------";

    cepgen::utils::Graph1D graph1d_bis("graph1d_bis", "cos(x)"), graph1d_ter("graph1d_ter", "cos(x)*x");
    for (const auto x : cepgen::Limits{-M_PI, M_PI}.generate(25)) {
      graph1d_bis.addPoint(x, cos(x));
      graph1d_ter.addPoint(x, cos(x) * x);
    }
    (void)plt->draw({&graph1d, &graph1d_bis, &graph1d_ter},
                    "multiplot1",
                    "a beautiful multiplot",
                    cepgen::utils::Drawer::Mode::grid);

    CG_LOG << "\n------- graph and hist -------";

    cepgen::utils::Hist1D hist1d_bis(graph1d.points().size(), {-M_PI, M_PI}, "hist1d_bis", "histogram");
    for (size_t i = 0; i < 10000; ++i)
      hist1d_bis.fill(gaussian1(gen));
    hist1d_bis.normalise(10.);
    (void)plt->draw({&graph1d, &hist1d_bis}, "multiplot2");

    cepgen::utils::Hist1D empty_hist(1, {0., 1.}, "empty histogram");
    (void)plt->draw(empty_hist);
    cepgen::utils::Graph1D empty_graph("empty graph");
    (void)plt->draw(empty_graph);
  }
  CG_TEST_SUMMARY;
}
