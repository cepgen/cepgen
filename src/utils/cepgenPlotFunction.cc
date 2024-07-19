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

#include "CepGen/Core/Exception.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Modules/FunctionalFactory.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Functional.h"
#include "CepGen/Utils/Graph.h"

using namespace std;

int main(int argc, char* argv[]) {
  string function, plotter;
  vector<string> functionals;
  int num_points;
  cepgen::Limits xrange, yrange;
  bool log, draw_grid, func_2d, ratio_plot;

  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("function,f", "function to parse", &function, "min(1.,exp(-x/10))")
      .addOptionalArgument("functional-eval,F", "functional evaluators", &functionals, vector<string>{})
      .addOptionalArgument("num-points,n", "number of points to consider", &num_points, 100)
      .addOptionalArgument("x-range,x", "horizontal axis range", &xrange, cepgen::Limits{-5., +5.})
      .addOptionalArgument("y-range,y", "vertical axis range", &yrange, cepgen::Limits{})
      .addOptionalArgument("draw-grid,g", "draw the x/y grid", &draw_grid, false)
      .addOptionalArgument("log,l", "logarithmic y-axis", &log, false)
      .addOptionalArgument("plotter,p", "type of plotter to user", &plotter, "")
      .addOptionalArgument("ratio,r", "draw the ratio plot", &ratio_plot, false)
      .addOptionalArgument("2d,t", "two-dimensional function", &func_2d, false)
      .parse();

  cepgen::initialise();

  if (functionals.empty())
    functionals = cepgen::FunctionalFactory::get().modules();
  CG_LOG << "Function to be plotted: " << function << ", functional evaluators: " << functionals << ".";

  // maps functional evaluators -> graphs
  map<string, cepgen::utils::Graph1D> m_gr_fb;
  map<string, cepgen::utils::Graph2D> m_gr2d_fb;
  for (const auto& func : functionals) {
    CG_LOG << "Building \"" << func << "\" functional.";
    try {
      vector<string> vars{"x"};
      if (func_2d)
        vars.emplace_back("y");
      auto test = cepgen::FunctionalFactory::get().build(
          func, cepgen::ParametersList().set<string>("expression", function).set<vector<string> >("variables", vars));
      if (func_2d) {
        if (!yrange.hasMin())
          yrange.min() = xrange.min();
        if (!yrange.hasMax())
          yrange.max() = xrange.max();
        for (int i = 0; i < num_points; ++i) {
          const double x = xrange.min() + (xrange.max() - xrange.min()) / (num_points - 1) * i;
          for (int j = 0; j < num_points; ++j) {
            const double y = yrange.min() + (yrange.max() - yrange.min()) / (num_points - 1) * j;
            m_gr2d_fb[func].addPoint(x, y, (*test)({x, y}));
          }
        }
      } else
        for (int i = 0; i < num_points; ++i) {
          const double x = xrange.min() + (xrange.max() - xrange.min()) / (num_points - 1) * i;
          m_gr_fb[func].addPoint(x, (*test)(x));
        }
    } catch (const cepgen::Exception&) {
      CG_WARNING("main") << "Exception encountered in \"" << func << "\" functional builder.";
      continue;
    }
  }

  if (!plotter.empty()) {
    auto plt = cepgen::DrawerFactory::get().build(plotter);
    cepgen::utils::Drawer::Mode dm;
    if (draw_grid)
      dm |= cepgen::utils::Drawer::Mode::grid;
    if (ratio_plot)
      dm |= cepgen::utils::Drawer::Mode::ratio;

    if (func_2d) {
      if (log)
        dm |= cepgen::utils::Drawer::Mode::logz;
      for (auto& gr_fb : m_gr2d_fb) {
        gr_fb.second.setName("graph2d_" + gr_fb.first);
        gr_fb.second.setTitle(gr_fb.first);
        plt->draw(gr_fb.second, dm);
      }
    } else {
      if (log)
        dm |= cepgen::utils::Drawer::Mode::logy;
      cepgen::utils::DrawableColl mg;
      for (auto& gr_fb : m_gr_fb) {
        gr_fb.second.setTitle(gr_fb.first);
        mg.emplace_back(&gr_fb.second);
      }
      plt->draw(mg, "comp_functionals", "", dm);
    }
  }

  return 0;
}
