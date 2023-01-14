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

#include <fstream>

#include "CepGen/CollinearFluxes/Parameterisation.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/Modes.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/String.h"

using namespace std;

int main(int argc, char* argv[]) {
  vector<string> cfluxes;
  int strfun_type, num_points;
  double mx;
  string ffmode, output_file, plotter;
  bool logx, logy, draw_grid, normalised;
  cepgen::Limits x_range, y_range, q2_range;

  cepgen::initialise();

  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument(
          "collflux,i", "collinear flux modelling(s)", &cfluxes, cepgen::CollinearFluxFactory::get().modules())
      .addOptionalArgument("formfac,f", "form factors modelling", &ffmode, "StandardDipole")
      .addOptionalArgument("mx,M", "diffractive mass (GeV/c^2)", &mx, 100.)
      .addOptionalArgument("sf,s", "structure functions modelling", &strfun_type, 301)
      .addOptionalArgument("xrange,x", "fractional loss range", &x_range, cepgen::Limits{0., 1.})
      .addOptionalArgument("yrange,y", "y range", &y_range)
      .addOptionalArgument("q2range,q", "parton virtuality range (GeV^2)", &q2_range, cepgen::Limits{0., 1000.})
      .addOptionalArgument("npoints,n", "number of x-points to scan", &num_points, 500)
      .addOptionalArgument("output,o", "output file name", &output_file, "collflux.scan.output.txt")
      .addOptionalArgument("plotter,p", "type of plotter to user", &plotter, "")
      .addOptionalArgument("logx", "logarithmic x-scale", &logx, false)
      .addOptionalArgument("logy,l", "logarithmic y-scale", &logy, false)
      .addOptionalArgument("draw-grid,g", "draw the x/y grid", &draw_grid, false)
      .addOptionalArgument("normalised", "plot xf(x) instead of f(x)", &normalised, false)
      .parse();

  ofstream out(output_file);
  if (logx && x_range.min() == 0.)
    x_range.min() = 1.e-3;
  const auto mx2 = mx * mx;

  out << "# coll. fluxes: " << cepgen::utils::merge(cfluxes, ",") << "\n"
      << "# struct. functions: " << strfun_type << "\n"
      << "# form factors: " << ffmode << "\n"
      << "# virtuality: " << q2_range << " GeV^2\n"
      << "# diffractive mass: " << mx << " GeV/c2\n"
      << "# fractional momentum loss: " << x_range;
  map<string, cepgen::utils::Graph1D> m_gr_fluxes;  // {collinear flux -> graph}
  map<string, unique_ptr<cepgen::collflux::Parameterisation> > coll_fluxes;
  for (const auto& cflux : cfluxes)
    coll_fluxes[cflux] = move(cepgen::CollinearFluxFactory::get().build(
        cflux,
        cepgen::ParametersList()
            .set("q2range", q2_range)
            .set("formFactors", ffmode)
            .set("structureFunctions",
                 cepgen::StructureFunctionsFactory::get().describeParameters(strfun_type).parameters())));

  for (const auto& x : x_range.generate(num_points, logx)) {
    vector<double> values;
    for (const auto& cf : coll_fluxes) {
      auto fx = (*cf.second)(x, mx2);
      if (isnan(fx))
        fx = -1.;
      values.emplace_back(fx);
      m_gr_fluxes[cf.first].addPoint(x, normalised ? x * fx : fx);
    }
    out << "\n" << x << "\t" << cepgen::utils::merge(values, "\t");
  }
  out.close();

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
    for (auto& cf_gr : m_gr_fluxes) {
      auto& gr = cf_gr.second;
      gr.xAxis().setLabel("x");
      gr.yAxis().setLabel(normalised ? "x dN/dx" : "dN/dx");
      if (y_range.valid())
        gr.yAxis().setRange(y_range);
      gr.setTitle(cepgen::utils::format("%s", cepgen::CollinearFluxFactory::get().describe(cf_gr.first).data()));
      coll.emplace_back(&gr);
    }
    plt->draw(coll, "comp_collfluxes", "", dm);
  }

  return 0;
}
