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
#include "CepGen/Integration/AnalyticIntegrator.h"
#include "CepGen/Modules/AnalyticIntegratorFactory.h"
#include "CepGen/Modules/CollinearFluxFactory.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/String.h"

using namespace std;

int main(int argc, char* argv[]) {
  vector<string> cfluxes;
  int num_points;
  double q2max, mxmin, mxmax;
  string integrator, output_file, plotter;
  bool logx, logy, draw_grid;
  vector<double> rescl, sqrts, accept;
  cepgen::Limits y_range;
  vector<cepgen::Limits> xi_ranges(1);

  cepgen::initialise();

  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("collflux,f",
                           "collinear flux modelling(s)",
                           &cfluxes,
                           cepgen::collflux::CollinearFluxFactory::get().modules())
      .addOptionalArgument("rescaling,r", "luminosity rescaling", &rescl, vector<double>{1.})
      .addOptionalArgument("integrator,i", "type of integration algorithm", &integrator, "gsl")
      .addOptionalArgument("q2max,q", "maximum Q^2", &q2max, 1000.)
      .addOptionalArgument("sqrts,s", "two-proton centre of mass energy (GeV)", &sqrts, vector<double>{13.e3})
      .addOptionalArgument("mxmin,m", "minimal two-photon mass", &mxmin, 0.)
      .addOptionalArgument("mxmax,M", "maximal two-photon mass", &mxmax, 100.)
      .addOptionalArgument("npoints,n", "number of x-points to scan", &num_points, 500)
      .addOptionalArgument("output,o", "output file name", &output_file, "collflux.int.scan.output.txt")
      .addOptionalArgument("plotter,p", "type of plotter to user", &plotter, "")
      .addOptionalArgument("logx", "logarithmic x-scale", &logx, false)
      .addOptionalArgument("logy,l", "logarithmic y-scale", &logy, false)
      .addOptionalArgument("xi-range,x", "acceptance range for proton momentum loss", &xi_ranges[0])
      .addOptionalArgument("accept", "pairs of min/max acceptance ranges for proton momentum loss", &accept)
      .addOptionalArgument("yrange,y", "y plot range", &y_range)
      .addOptionalArgument("draw-grid,g", "draw the x/y grid", &draw_grid, false)
      .parse();

  ofstream out(output_file);
  if (logx && mxmin == 0.)
    mxmin = 1.e-3;
  if (sqrts.size() != cfluxes.size())
    sqrts = vector<double>(cfluxes.size(), sqrts.at(0));
  if (rescl.size() != cfluxes.size())
    rescl = vector<double>(cfluxes.size(), rescl.at(0));
  if (!accept.empty()) {
    xi_ranges.clear();
    if (accept.size() % 2 != 0)
      throw CG_FATAL("main")
          << "Invalid acceptance(s) list specified! Supported format is an array of (min1, max2, min2, max2, ...)";
    for (size_t i = 0; i < accept.size() / 2; ++i) {
      cepgen::Limits rng{accept.at(2 * i), accept.at(2 * i + 1)};
      if (rng == cepgen::Limits{0., 1.})
        rng = cepgen::Limits();
      xi_ranges.emplace_back(rng);
    }
    if (!xi_ranges.empty())
      CG_LOG << "x (xi) acceptance cuts defined: " << xi_ranges << ".";
  }

  out << "# coll. fluxes: " << cepgen::utils::merge(cfluxes, ",") << "\n"
      << "# two-photon mass range: " << cepgen::Limits(mxmin, mxmax);
  map<string, vector<cepgen::utils::Graph1D> > m_gr_fluxes;  // {collinear flux -> graph}
  vector<double> mxvals;
  for (int j = 0; j < num_points; ++j)
    mxvals.emplace_back((!logx) ? mxmin + j * (mxmax - mxmin) / (num_points - 1)
                                : pow(10, log10(mxmin) + j * (log10(mxmax) - log10(mxmin)) / (num_points - 1)));
  vector<vector<double> > values(num_points);
  auto integr = cepgen::AnalyticIntegratorFactory::get().build(
      cepgen::ParametersList().setName<string>(integrator).set<int>("mode", 0).set<int>("nodes", 2000));
  for (size_t i = 0; i < cfluxes.size(); ++i) {
    const auto& cflux = cfluxes.at(i);
    const auto s = sqrts.at(i) * sqrts.at(i);
    auto coll_flux = cepgen::collflux::CollinearFluxFactory::get().build(cflux);
    for (const auto& xi_range : xi_ranges) {
      ostringstream oss;
      oss << cflux;
      if (xi_range.valid())
        oss << " (" << xi_range.min() << " < \\xi < " << xi_range.max() << ")";
      m_gr_fluxes[cflux].emplace_back();
      m_gr_fluxes[cflux].back().setTitle(oss.str());
    }

    for (int j = 0; j < num_points; ++j) {
      const auto& mx = mxvals.at(j);
      for (size_t k = 0; k < xi_ranges.size(); ++k) {
        const auto& xi_range = xi_ranges.at(k);
        auto lumi_wgg = integr->eval(
            [&xi_range, &mx, &s, &coll_flux](double x) {
              if (xi_range.valid() && (!xi_range.contains(x) || !xi_range.contains(mx * mx / x / s)))
                return 0.;
              return 2. * mx / x / s * (*coll_flux)(x) * (*coll_flux)(mx * mx / x / s);
            },
            cepgen::Limits(mx * mx / s, 1.));
        lumi_wgg *= rescl.at(i);
        values.at(j).emplace_back(lumi_wgg);
        m_gr_fluxes[cflux][k].addPoint(mx, lumi_wgg);
      }
    }
  }
  for (int i = 0; i < num_points; ++i)
    out << "\n" << mxvals.at(i) << "\t" << cepgen::utils::merge(values.at(i), "\t");
  out.close();

  if (!plotter.empty()) {
    auto plt = cepgen::utils::DrawerFactory::get().build(plotter);
    cepgen::utils::Drawer::Mode dm;
    if (logx)
      dm |= cepgen::utils::Drawer::Mode::logx;
    if (logy)
      dm |= cepgen::utils::Drawer::Mode::logy;
    if (draw_grid)
      dm |= cepgen::utils::Drawer::Mode::grid;
    cepgen::utils::DrawableColl coll;
    size_t i = 0;
    for (auto& cf_gr : m_gr_fluxes) {
      for (auto& gr_xi : cf_gr.second) {
        string units;
        if (rescl.at(i) != 1.)
          units = "cm$^{-2}$ s$^{-1}$ ";
        gr_xi.xAxis().setLabel("$\\omega_{\\gamma\\gamma}$ (GeV)");
        gr_xi.yAxis().setLabel("d$L_{\\gamma\\gamma}$/d$\\omega_{\\gamma\\gamma}$ (" + units + "GeV$^{-1}$)");
        if (y_range.valid())
          gr_xi.yAxis().setRange(y_range);
        //gr_xi.setTitle(
        //    cepgen::utils::format("%s", cepgen::collflux::CollinearFluxFactory::get().describe(cf_gr.first).data()));
        coll.emplace_back(&gr_xi);
      }
      ++i;
    }
    plt->draw(coll, "comp_photonlumi", "", dm);
  }
  return 0;
}
