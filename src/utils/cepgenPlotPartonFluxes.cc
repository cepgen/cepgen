/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2023  Laurent Forthomme
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

#include "CepGen/Core/Exception.h"
#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Modules/FormFactorsFactory.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/Beam.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/String.h"

using namespace std;

int main(int argc, char* argv[]) {
  vector<string> fluxes_names;
  int num_points;
  double kt2, mx;
  string output_file, plotter;
  bool logx, logy, draw_grid, normalised;
  cepgen::Limits x_range, y_range, q2_range;

  cepgen::initialise();

  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument(
          "fluxes", "parton fluxe modellings", &fluxes_names, cepgen::PartonFluxFactory::get().modules())
      .addOptionalArgument("mx,M", "diffractive mass (GeV)", &mx, 1.5)
      .addOptionalArgument("xrange,x", "fractional loss range", &x_range, cepgen::Limits{0., 1.})
      .addOptionalArgument("yrange,y", "y range", &y_range)
      .addOptionalArgument("q2range,q", "parton virtuality range (GeV^2)", &q2_range, cepgen::Limits{0., 1000.})
      .addOptionalArgument("kt2,k", "parton transverse virtuality (GeV^2)", &kt2, 100.)
      .addOptionalArgument("npoints,n", "number of x-points to scan", &num_points, 100)
      .addOptionalArgument("output,o", "output file name", &output_file, "flux.scan.output.txt")
      .addOptionalArgument("plotter,p", "type of plotter to user", &plotter, "")
      .addOptionalArgument("logx", "logarithmic x-scale", &logx, false)
      .addOptionalArgument("logy,l", "logarithmic y-scale", &logy, false)
      .addOptionalArgument("draw-grid,g", "draw the x/y grid", &draw_grid, false)
      .addOptionalArgument("normalised", "plot xf(x) instead of f(x)", &normalised, false)
      .parse();

  const double mx2 = mx * mx;
  if (logx && x_range.min() == 0.)
    x_range.min() = 1.e-3;
  if (x_range.max() == 1.)
    x_range.max() -= 1.e-15;

  ofstream out(output_file);
  out << "# parton fluxes: " << cepgen::utils::merge(fluxes_names, ",") << "\n"
      << "# virtuality: " << kt2 << " GeV^2\n"
      << "# diffractive mass: " << mx << " GeV/c2\n"
      << "# fractional momentum loss: " << x_range;

  vector<std::unique_ptr<cepgen::PartonFlux> > fluxes;
  vector<cepgen::utils::Graph1D> graph_flux;
  for (const auto& flux : fluxes_names) {
    fluxes.emplace_back(cepgen::PartonFluxFactory::get().build(flux));
    graph_flux.emplace_back(flux, cepgen::PartonFluxFactory::get().describe(flux));
  }
  for (const auto& x : x_range.generate(num_points)) {
    out << "\n" << x;
    double f = 0.;
    for (const auto& flux_vs_name : fluxes) {
      if (flux_vs_name.second->ktFactorised())
        f = (*flux_vs_name.second)(x, kt2, mx2);
      else
        f = (*flux_vs_name.second)(x, mx2);
      if (normalised)
        f *= x;
      out << "\t" << f;
      graph_flux.at(flux_vs_name.first).addPoint(x, f);
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
    cepgen::utils::DrawableColl coll;

    for (auto& gr : graph_flux) {
      gr.xAxis().setLabel("$\\xi$");
      gr.yAxis().setLabel(normalised ? "$\\xi\\varphi_{T}(\\xi, k_{T}^{2})" : "$\\varphi_{T}(\\xi, k_{T}^{2})$");
      if (y_range.valid())
        gr.yAxis().setRange(y_range);
      gr.setTitle(cepgen::utils::format("%s", cepgen::PartonFluxFactory::get().describe(gr.first).data()));
      coll.emplace_back(&gr);
    }
    plt->draw(coll, "comp_partonflux", cepgen::utils::format("$k_{T}^{2}$ = %g GeV$^{2}$", kt2), dm);
  }

  return 0;
}
