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

#include "CepGen/Generator.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/Histogram.h"

using namespace std;

int main(int argc, char* argv[]) {
  string plotter;

  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("plotter,p", "type of plotter to user", &plotter, "text")
      .parse();
  cepgen::initialise();

  auto plt = cepgen::utils::DrawerFactory::get().build(plotter);

  cout << "---------- 1D graph ----------" << endl;

  // test 1D graph
  cepgen::utils::Graph1D graph1d("graph1d", "sin(x)");
  for (double x = -M_PI; x <= M_PI; x += 0.25)
    graph1d.addPoint(x, sin(x));
  plt->draw(graph1d);

  cout << endl << "---------- 2D graph ----------" << endl;

  // test 2D graph
  cepgen::utils::Graph2D graph2d("graph2d");
  for (double x = -5.; x < 5.; x += 0.5)
    for (double y = -5.; y < 5.; y += 0.2)
      graph2d.addPoint(x, y, (sin(x) / x) * (sin(y) / y));
  plt->draw(graph2d);

  default_random_engine gen;

  cout << endl << "-------- 1D histogram --------" << endl;

  // test 1D histogram
  cepgen::utils::Hist1D hist1d(20, {-5., 5.}, "hist1d");
  cauchy_distribution<double> bw(0., 1.);
  for (size_t i = 0; i < 10000; ++i)
    hist1d.fill(bw(gen));
  hist1d.xAxis().setLabel("Random variable");
  hist1d.yAxis().setLabel("Occurrences");
  plt->draw(hist1d);

  cout << endl << "-------- 2D histogram --------" << endl;

  // test 2d histogram
  cepgen::utils::Hist2D hist2d(20, {-5., 5.}, 50, {-5., 5.}, "hist2d");
  normal_distribution<double> gaus1(0., 1.), gaus2(0., 1.);
  for (size_t i = 0; i < 1000; ++i)
    for (size_t j = 0; j < 1000; ++j)
      hist2d.fill(gaus1(gen), gaus2(gen));
  plt->draw(hist2d, cepgen::utils::Drawer::Mode::logz);

  cout << endl << "--------- multiplots ---------" << endl;

  cepgen::utils::Graph1D graph1d_bis("graph1d_bis", "cos(x)"), graph1d_ter("graph1d_ter", "cos(x)*x");
  for (double x = -M_PI; x <= M_PI; x += 0.25) {
    graph1d_bis.addPoint(x, cos(x));
    graph1d_ter.addPoint(x, cos(x) * x);
  }
  plt->draw({&graph1d, &graph1d_bis, &graph1d_ter}, "multiplot1", "a beautiful multiplot");

  cout << endl << "------- graph and hist -------" << endl;

  cepgen::utils::Hist1D hist1d_bis(graph1d.points().size(), {-M_PI, M_PI}, "hist1d_bis");
  for (size_t i = 0; i < 10000; ++i)
    hist1d_bis.fill(gaus1(gen));
  hist1d_bis.normalise(10.);
  plt->draw({&graph1d, &hist1d_bis}, "multiplot2");

  return 0;
}
