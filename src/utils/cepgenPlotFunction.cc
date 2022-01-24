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
  int num_points;
  double min_x, max_x;
  bool logy, draw_grid;

  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("function,f", "function to parse", &function, "min(1.,exp(-x/10))")
      .addOptionalArgument("num-points,n", "number of points to consider", &num_points, 100)
      .addOptionalArgument("min-x,m", "minimal range", &min_x, -5.)
      .addOptionalArgument("max-x,M", "maximal range", &max_x, +5.)
      .addOptionalArgument("draw-grid,g", "draw the x/y grid", &draw_grid, false)
      .addOptionalArgument("logy,l", "logarithmic y-axis", &logy, false)
      .addOptionalArgument("plotter,p", "type of plotter to user", &plotter, "")
      .parse();

  cepgen::initialise();

  CG_LOG << "Function to be plotted: " << function;

  map<string, cepgen::utils::Graph1D> m_gr_fb;
  for (const auto& func : cepgen::utils::FunctionalFactory::get().modules()) {
    CG_LOG << "Building \"" << func << "\" functional.";
    try {
      auto test =
          cepgen::utils::FunctionalFactory::get().build(func,
                                                        cepgen::ParametersList()
                                                            .set<std::string>("expression", function)
                                                            .set<std::vector<std::string> >("variables", {"x"}));
      for (unsigned short i = 0; i < num_points; ++i) {
        const double x = min_x + (max_x - min_x) / (num_points - 1) * i;
        m_gr_fb[func].addPoint(x, (*test)(x));
      }
    } catch (const cepgen::Exception&) {
      CG_WARNING("main") << "Exception encountered in \"" << func << "\" functional builder.";
      continue;
    }
  }

  if (!plotter.empty()) {
    auto plt = cepgen::utils::DrawerFactory::get().build(plotter);
    cepgen::utils::Drawer::Mode dm;
    if (logy)
      dm |= cepgen::utils::Drawer::Mode::logy;
    if (draw_grid)
      dm |= cepgen::utils::Drawer::Mode::grid;

    cepgen::utils::DrawableColl mg;
    for (auto& gr_fb : m_gr_fb) {
      gr_fb.second.setTitle(gr_fb.first);
      mg.emplace_back(&gr_fb.second);
    }
    plt->draw(mg, "comp_functionals", "", dm);
  }

  return 0;
}
