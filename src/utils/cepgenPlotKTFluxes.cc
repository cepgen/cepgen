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

#include "CepGen/Generator.h"
#include "CepGen/KTFluxes/KTFlux.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/Message.h"
#include "CepGen/Utils/String.h"

using namespace std;
using namespace std::string_literals;

int main(int argc, char* argv[]) {
  vector<string> fluxes_names;
  int num_points;
  double kt2, mx, q2;
  string output_file, plotter;
  bool logx, logy, draw_grid, normalised, ratio_plot;
  cepgen::Limits x_range, y_range;

  cepgen::initialise();

  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument(
          "fluxes,f", "parton fluxes modellings", &fluxes_names, cepgen::KTFluxFactory::get().modules())
      .addOptionalArgument("mx,M", "diffractive mass (GeV)", &mx, 1.5)
      .addOptionalArgument("q2,q", "parton virtuality (GeV^2)", &q2, -1.)
      .addOptionalArgument("xrange,x", "fractional loss range", &x_range, cepgen::Limits{0., 1.})
      .addOptionalArgument("yrange,y", "y range", &y_range)
      .addOptionalArgument("kt2,k", "parton transverse virtuality (GeV^2)", &kt2, 10.)
      .addOptionalArgument("npoints,n", "number of x-points to scan", &num_points, 100)
      .addOptionalArgument("output,o", "output file name", &output_file, "flux.scan.output.txt")
      .addOptionalArgument("plotter,p", "type of plotter to user", &plotter, "")
      .addOptionalArgument("logx", "logarithmic x-scale", &logx, false)
      .addOptionalArgument("logy,l", "logarithmic y-scale", &logy, false)
      .addOptionalArgument("draw-grid,g", "draw the x/y grid", &draw_grid, false)
      .addOptionalArgument("ratio,r", "draw the ratio plot", &ratio_plot, false)
      .addOptionalArgument("normalised", "plot xf(x) instead of f(x)", &normalised, false)
      .parse();

  const bool plot_vs_q2 = (q2 >= 0.);
  if (logx && x_range.min() == 0.)
    x_range.min() = 1.e-3;
  if (x_range.max() == 1.)
    x_range.max() -= 1.e-15;

  vector<std::unique_ptr<cepgen::KTFlux> > fluxes;
  vector<cepgen::utils::Graph1D> graph_flux;
  for (const auto& flux : fluxes_names) {
    fluxes.emplace_back(cepgen::KTFluxFactory::get().build(flux));
    graph_flux.emplace_back(flux, cepgen::KTFluxFactory::get().describe(flux));
  }

  ofstream out(output_file);
  out << "# parton fluxes: " << cepgen::utils::merge(fluxes_names, ";") << "\n";
  out << "# transverse virtuality: " << kt2 << " GeV^2\n";
  if (plot_vs_q2)
    out << "# virtuality: " << q2 << " GeV^2\n";
  else
    out << "# diffractive mass: " << mx << " GeV/c^2\n";
  out << "# fractional momentum loss: " << x_range;

  for (const auto& x : x_range.generate(num_points)) {
    out << "\n" << x;
    for (size_t i = 0; i < fluxes.size(); ++i) {
      auto flux = plot_vs_q2 ? fluxes.at(i)->fluxQ2(x, kt2, q2) : fluxes.at(i)->fluxMX2(x, kt2, mx * mx);
      flux *= (normalised ? x : 1.);
      out << "\t" << flux;
      graph_flux.at(i).addPoint(x, flux);
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
      gr.xAxis().setLabel("$\\xi$");
      gr.yAxis().setLabel("$"s + (normalised ? "\\xi\\varphi" : "\\varphi") + "(\\xi, " +
                          (plot_vs_q2 ? "Q^{2}" : "M_{X}") + ", k_{T}^{2})$");
      if (y_range.valid())
        gr.yAxis().setRange(y_range);
      coll.emplace_back(&gr);
    }
    plt->draw(coll,
              "comp_partonflux",
              (plot_vs_q2 ? cepgen::utils::format("$Q^{2}$ = %g GeV$^{2}$", q2)
                          : cepgen::utils::format("$M_{X}$ = %g GeV", mx)) +
                  ", " + cepgen::utils::format("$k_{T}^{2}$ = %g GeV$^{2}$", kt2),
              dm);
  }

  return 0;
}
