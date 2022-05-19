/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022  Laurent Forthomme
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
#include "CepGen/Generator.h"
#include "CepGen/Modules/CollinearFluxFactory.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/GSLIntegrator.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/String.h"

using namespace std;

int main(int argc, char* argv[]) {
  vector<int> cfluxes;
  int num_points;
  double q2max, sqrts, mxmin, mxmax;
  string output_file, plotter;
  bool logx, logy, draw_grid;
  vector<double> xi_range_vec;
  cepgen::Limits y_range, xi_range;

  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("collflux,i", "collinear flux modelling(s)", &cfluxes, vector<int>{1})
      .addOptionalArgument("q2max,q", "maximum Q^2", &q2max, 1000.)
      .addOptionalArgument("sqrts,s", "two-proton centre of mass energy (GeV)", &sqrts, 13.e3)
      .addOptionalArgument("mxmin,m", "minimal two-photon mass", &mxmin, 0.)
      .addOptionalArgument("mxmax,M", "maximal two-photon mass", &mxmax, 100.)
      .addOptionalArgument("npoints,n", "number of x-points to scan", &num_points, 500)
      .addOptionalArgument("output,o", "output file name", &output_file, "collflux.int.scan.output.txt")
      .addOptionalArgument("plotter,p", "type of plotter to user", &plotter, "")
      .addOptionalArgument("logx", "logarithmic x-scale", &logx, false)
      .addOptionalArgument("logy,l", "logarithmic y-scale", &logy, false)
      .addOptionalArgument("xi-range,x", "acceptance range for proton momentum loss", &xi_range)
      .addOptionalArgument("yrange,y", "y plot range", &y_range)
      .addOptionalArgument("draw-grid,g", "draw the x/y grid", &draw_grid, false)
      .parse();

  cepgen::initialise();

  ofstream out(output_file);
  if (logx && mxmin == 0.)
    mxmin = 1.e-3;

  out << "# coll. fluxes: " << cepgen::utils::merge(cfluxes, ",") << "\n"
      << "# two-photon mass range: " << cepgen::Limits(mxmin, mxmax);
  map<int, cepgen::utils::Graph1D> m_gr_fluxes;  // {collinear flux -> graph}
  vector<double> mxvals;
  for (int j = 0; j < num_points; ++j)
    mxvals.emplace_back((!logx) ? mxmin + j * (mxmax - mxmin) / (num_points - 1)
                                : pow(10, log10(mxmin) + j * (log10(mxmax) - log10(mxmin)) / (num_points - 1)));
  vector<vector<double> > values(num_points);
  auto integr = cepgen::utils::GSLIntegrator();
  const auto s = sqrts * sqrts;
  for (const auto& cflux : cfluxes) {
    auto coll_flux = cepgen::collflux::CollinearFluxFactory::get().build(cflux);
    ostringstream oss;
    oss << cflux;
    m_gr_fluxes[cflux].setTitle(oss.str());

    for (int j = 0; j < num_points; ++j) {
      const auto& mx = mxvals.at(j);
      auto lumi_wgg = integr.eval(
          [&xi_range, &mx, &s, &coll_flux](double x) {
            if (xi_range.valid() && (!xi_range.contains(x) || !xi_range.contains(mx * mx / x / s)))
              return 0.;
            return 2. * mx / x / s * (*coll_flux)(x) * (*coll_flux)(mx * mx / x / s);
          },
          mx * mx / s,
          1.);
      values.at(j).emplace_back(lumi_wgg);
      m_gr_fluxes[cflux].addPoint(mx, lumi_wgg);
    }
  }
  for (int i = 0; i < num_points; ++i)
    out << "\n" << mxvals.at(i) << "\t" << cepgen::utils::merge(values.at(i), "\t");
  out.close();

  if (!plotter.empty()) {
    ostringstream oss;
    auto plt = cepgen::utils::DrawerFactory::get().build(plotter);
    cepgen::utils::Drawer::Mode dm;
    if (logx)
      dm |= cepgen::utils::Drawer::Mode::logx;
    if (logy)
      dm |= cepgen::utils::Drawer::Mode::logy;
    if (draw_grid)
      dm |= cepgen::utils::Drawer::Mode::grid;
    cepgen::utils::DrawableColl coll;
    for (auto& cf_gr : m_gr_fluxes) {
      cf_gr.second.xAxis().setLabel("$\\omega_{\\gamma\\gamma}$ (GeV$^{-1}$)");
      cf_gr.second.yAxis().setLabel("d$L_{\\gamma\\gamma}$/d$\\omega_{\\gamma\\gamma}$ (GeV$^{-1}$)");
      if (y_range.valid())
        cf_gr.second.yAxis().setRange(y_range);
      cf_gr.second.setTitle(
          cepgen::utils::format("%s", cepgen::collflux::CollinearFluxFactory::get().describe(cf_gr.first).data()));
      coll.emplace_back(&cf_gr.second);
    }
    plt->draw(coll, "comp_photonlumi", oss.str(), dm);
  }

  return 0;
}
