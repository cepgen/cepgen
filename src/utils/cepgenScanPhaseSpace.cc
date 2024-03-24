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
#include "CepGen/Core/RunParameters.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/String.h"

using namespace std;

int main(int argc, char* argv[]) {
  string input_card, plotter;
  int npoints;
  vector<int> dim;
  double def;
  bool draw_grid, log;

  cepgen::ArgumentsParser(argc, argv)
      .addArgument("input,i", "input card", &input_card)
      .addOptionalArgument("default,D", "default value for non-varying coordinates", &def, 0.5)
      .addOptionalArgument("dim,s", "dimensions to probe", &dim, vector<int>{0, 1})
      .addOptionalArgument("num-points,n", "number of points to probe", &npoints, 100)
      .addOptionalArgument("draw-grid,g", "draw the x/y grid", &draw_grid, false)
      .addOptionalArgument("log,l", "logarithmic axis", &log, false)
      .addOptionalArgument("plotter,p", "type of plotter to user", &plotter, "")
      .parse();

  cepgen::utils::Graph1D gr_scan_1d("test_scan");
  cepgen::utils::Graph2D gr_scan_2d("test_scan");
  if (dim.size() > 3)
    throw CG_FATAL("main") << "Number of dimensions to probe (" << dim.size() << ") is too high";

  cepgen::Generator gen;
  gen.parseRunParameters(input_card);
  CG_LOG << gen.runParameters();
  const size_t ndim = gen.runParameters().process().ndim();

  vector<double> coord(ndim, def);

  for (int i = 0; i < npoints; ++i) {
    const double x = i * 1. / npoints;
    switch (dim.size()) {
      case 0:
        gr_scan_1d.addPoint(x, gen.computePoint(vector<double>(ndim, x)));
        break;
      case 1:
        coord[dim.at(0)] = x;
        gr_scan_1d.addPoint(x, gen.computePoint(coord));
        break;
      case 2:
        coord[dim.at(0)] = x;
        for (int j = 0; j < npoints; ++j) {
          const double y = j * 1. / npoints;
          coord[dim.at(1)] = y;
          gr_scan_2d.addPoint(x, y, gen.computePoint(coord));
        }
        break;
    }
  }
  if (!plotter.empty()) {
    auto plt = cepgen::DrawerFactory::get().build(plotter);
    cepgen::utils::Drawer::Mode dm;
    if (draw_grid)
      dm |= cepgen::utils::Drawer::Mode::grid;
    switch (dim.size()) {
      case 0:
      case 1: {
        if (log)
          dm |= cepgen::utils::Drawer::Mode::logy;
        string xlabel;
        if (dim.empty())
          xlabel = cepgen::utils::format("x_{i = 0, ..., %ld}", ndim - 1);
        else
          xlabel = cepgen::utils::format("x_{%d}", dim.at(0));
        gr_scan_1d.setTitle(cepgen::utils::format("%s variation, all others x_{i} at %g", xlabel.c_str(), def));
        gr_scan_1d.xAxis().setLabel(xlabel);
        gr_scan_2d.yAxis().setLabel(cepgen::utils::format("d^{%ld}#sigma/d#bf{x}^{%ld}", ndim, ndim));
        plt->draw(gr_scan_1d, dm);
      } break;
      case 2: {
        if (log)
          dm |= cepgen::utils::Drawer::Mode::logz;
        string xlabel = cepgen::utils::format("x_{%d}", dim.at(0)), ylabel = cepgen::utils::format("x_{%d}", dim.at(1));
        gr_scan_2d.setTitle(
            cepgen::utils::format("(%s, %s) variation, all others x_{i} at %g", xlabel.c_str(), ylabel.c_str(), def));
        gr_scan_2d.xAxis().setLabel(xlabel);
        gr_scan_2d.yAxis().setLabel(ylabel);
        gr_scan_2d.zAxis().setLabel(cepgen::utils::format("d^{%ld}#sigma/d#bf{x}^{%ld}", ndim, ndim));
        plt->draw(gr_scan_2d, dm);
      }
    }
  }

  return 0;
}
