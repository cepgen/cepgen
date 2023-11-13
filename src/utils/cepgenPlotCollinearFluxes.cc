/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023  Laurent Forthomme
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

#include <fstream>

#include "CepGen/CollinearFluxes/CollinearFlux.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/String.h"

using namespace std;
using namespace std::string_literals;

int main(int argc, char* argv[]) {
  vector<string> fluxes_names;
  int num_points;
  double mx, q2;
  string output_file, plotter;
  bool logx, logy, draw_grid, normalised, ratio_plot;
  cepgen::Limits x_range, y_range;

  cepgen::initialise();

  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument(
          "fluxes,f", "parton fluxes modellings", &fluxes_names, cepgen::CollinearFluxFactory::get().modules())
      .addOptionalArgument("mx,M", "diffractive mass (GeV)", &mx, 1.5)
      .addOptionalArgument("q2,q", "parton virtuality (GeV^2)", &q2, -1.)
      .addOptionalArgument("xrange,x", "fractional loss range", &x_range, cepgen::Limits{0., 1.})
      .addOptionalArgument("yrange,y", "y range", &y_range)
      .addOptionalArgument("npoints,n", "number of x-points to scan", &num_points, 100)
      .addOptionalArgument("output,o", "output file name", &output_file, "flux.scan.output.txt")
      .addOptionalArgument("plotter,p", "type of plotter to user", &plotter, "")
      .addOptionalArgument("logx", "logarithmic x-scale", &logx, false)
      .addOptionalArgument("logy,l", "logarithmic y-scale", &logy, false)
      .addOptionalArgument("draw-grid,g", "draw the x/y grid", &draw_grid, false)
      .addOptionalArgument("ratio,r", "draw the ratio plot", &ratio_plot, false)
      .addOptionalArgument("normalised", "plot xf(x) instead of f(x)", &normalised, false)
      .parse();

  const bool plot_vs_q2 = (q2 > 0.);
  const double mx2 = mx * mx;
  if (logx && x_range.min() == 0.)
    x_range.min() = 1.e-3;
  if (x_range.max() == 1.)
    x_range.max() -= 1.e-15;

  vector<std::unique_ptr<cepgen::CollinearFlux> > fluxes;
  vector<cepgen::utils::Graph1D> graph_flux;
  for (const auto& flux : fluxes_names) {
    fluxes.emplace_back(cepgen::CollinearFluxFactory::get().build(flux));
    graph_flux.emplace_back(flux, cepgen::CollinearFluxFactory::get().describe(flux));
  }

  ofstream out(output_file);
  out << "# parton fluxes: " << cepgen::utils::merge(fluxes_names, ";") << "\n";
  out << "# fractional momentum loss: " << x_range;

  for (const auto& x : x_range.generate(num_points)) {
    out << "\n" << x;
    for (size_t j = 0; j < fluxes.size(); ++j) {
      auto flux = plot_vs_q2 ? fluxes.at(j)->fluxQ2(x, q2) : fluxes.at(j)->fluxMX2(x, mx2);
      flux *= (normalised ? x : 1.);
      out << "\t" << flux;
      graph_flux.at(j).addPoint(x, flux);
    }
  }
  out.close();
  CG_LOG << "Scan written in \"" << output_file << "\".";

  if (!plotter.empty()) {
    auto plt = cepgen::DrawerFactory::get().build(plotter);
    cepgen::utils::Drawer::Mode dm;
    if (logx)
      dm |= cepgen::utils::Drawer::Mode::logx;
    if (logy)
      dm |= cepgen::utils::Drawer::Mode::logy;
    if (draw_grid)
      dm |= cepgen::utils::Drawer::Mode::grid;
    if (ratio_plot)
      dm |= cepgen::utils::Drawer::Mode::ratio;
    cepgen::utils::DrawableColl coll;

    for (auto& gr : graph_flux) {
      gr.xAxis().setLabel("$x$");
      gr.yAxis().setLabel("$"s + (normalised ? "xf" : "f") + "(x," + (plot_vs_q2 ? "Q^{2}" : "M_{X}") + ")$");
      if (y_range.valid())
        gr.yAxis().setRange(y_range);
      coll.emplace_back(&gr);
    }
    plt->draw(coll,
              "comp_partonflux",
              plot_vs_q2 ? cepgen::utils::format("$Q^{2}$ = %g GeV$^{2}$", q2)
                         : cepgen::utils::format("$M_{X}$ = %g GeV/c$^{2}$", mx),
              dm);
  }

  return 0;
}
