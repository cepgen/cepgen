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
#include "CepGen/Modules/CollinearFluxFactory.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/Modes.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/String.h"

using namespace std;
using namespace cepgen;

int main(int argc, char* argv[]) {
  vector<int> cfluxes, modes;
  int strfun_type, num_points;
  double mx, xmin, xmax;
  string ffmode, output_file, plotter;
  bool logx, logy, draw_grid;

  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("collflux,i", "collinear flux modelling(s)", &cfluxes, vector<int>{1})
      .addOptionalArgument("formfac,f", "form factors modelling", &ffmode, "StandardDipole")
      .addOptionalArgument("modes,t", "beam modelling(s)", &modes, vector<int>{(int)cepgen::Beam::Mode::ProtonElastic})
      .addOptionalArgument("mx,M", "diffractive mass (GeV/c^2)", &mx, 100.)
      .addOptionalArgument("sf,s", "structure functions modelling", &strfun_type, 301)
      .addOptionalArgument("xmin,x", "minimal fractional loss", &xmin, 0.)
      .addOptionalArgument("xmax,X", "maximal fractional loss", &xmax, 1.)
      .addOptionalArgument("npoints,n", "number of x-points to scan", &num_points, 500)
      .addOptionalArgument("output,o", "output file name", &output_file, "collflux.scan.output.txt")
      .addOptionalArgument("plotter,p", "type of plotter to user", &plotter, "")
      .addOptionalArgument("logx", "logarithmic x-scale", &logx, false)
      .addOptionalArgument("logy,l", "logarithmic y-scale", &logy, false)
      .addOptionalArgument("draw-grid,g", "draw the x/y grid", &draw_grid, false)
      .parse();

  cepgen::initialise();

  ofstream out(output_file);
  if (logx && xmin == 0.)
    xmin = 1.e-3;

  out << "# coll. fluxes: " << cepgen::utils::merge(cfluxes, ",") << "\n"
      << "# struct. functions: " << strfun_type << "\n"
      << "# form factors: " << ffmode << "\n"
      << "# diffractive mass: " << mx << " GeV/c2\n"
      << "# fractional momentum loss: " << cepgen::Limits(xmin, xmax) << "\n"
      << "# fluxes modes: " << cepgen::utils::merge(modes, ",") << "\n";
  map<int, vector<cepgen::utils::Graph1D> > m_v_gr_fluxes;  // {collinear flux -> {mode, mode, ...}}
  for (const auto& cflux : cfluxes) {
    vector<unique_ptr<cepgen::collflux::Parameterisation> > coll_fluxes;
    for (const auto& mode : modes) {
      cepgen::ParametersList flux_params;
      flux_params.set("q2range", cepgen::Limits{0., 1000.})
          .set("formFactors", ffmode)
          .set("structureFunctions",
               cepgen::strfun::StructureFunctionsFactory::get().describeParameters(strfun_type).parameters());
      switch ((cepgen::Beam::Mode)mode) {
        case cepgen::Beam::Mode::ProtonElastic:
          flux_params.setAs<int, cepgen::Beam::KTFlux>("ktFlux", cepgen::Beam::KTFlux::P_Photon_Elastic);
          break;
        case cepgen::Beam::Mode::ProtonInelastic:
          flux_params.setAs<int, cepgen::Beam::KTFlux>("ktFlux", cepgen::Beam::KTFlux::P_Photon_Inelastic);
          break;
        default:
          throw CG_FATAL("main") << "Invalid beam mode: " << mode << "!";
      }
      coll_fluxes.emplace_back(cepgen::collflux::CollinearFluxFactory::get().build(cflux, flux_params));
      m_v_gr_fluxes[cflux].emplace_back();
      ostringstream oss;
      oss << (cepgen::Beam::Mode)mode;
      m_v_gr_fluxes[cflux].back().setTitle(oss.str());
    }

    for (int i = 0; i < num_points; ++i) {
      const double x = (!logx) ? xmin + i * (xmax - xmin) / (num_points - 1)
                               : pow(10, log10(xmin) + i * (log10(xmax) - log10(xmin)) / (num_points - 1));
      out << x;
      for (size_t j = 0; j < coll_fluxes.size(); ++j) {
        auto fx = (*coll_fluxes.at(j))(x, mx);
        if (isnan(fx))
          fx = -1.;
        out << "\t" << fx;
        m_v_gr_fluxes[cflux].at(j).addPoint(x, fx);
      }
      out << "\n";
    }
  }
  out.close();

  if (!plotter.empty()) {
    ostringstream oss;
    oss << "M_{X} = " << mx << " GeV/c^{2} (" << ffmode << ", " << (cepgen::strfun::Type)strfun_type << ")";
    auto plt = cepgen::utils::DrawerFactory::get().build(plotter);
    cepgen::utils::Drawer::Mode dm;
    if (logx)
      dm |= cepgen::utils::Drawer::Mode::logx;
    if (logy)
      dm |= cepgen::utils::Drawer::Mode::logy;
    if (draw_grid)
      dm |= cepgen::utils::Drawer::Mode::grid;
    cepgen::utils::DrawableColl coll;
    for (auto& cf_gr : m_v_gr_fluxes)
      for (auto& gr : cf_gr.second) {
        gr.xAxis().setLabel("x");
        gr.yAxis().setLabel("dN/dx");
        if (cf_gr.second.size() > 1)
          gr.setTitle(cepgen::utils::format("%s - %s",
                                            cepgen::collflux::CollinearFluxFactory::get().describe(cf_gr.first).data(),
                                            gr.title().data()));
        else
          gr.setTitle(
              cepgen::utils::format("%s", cepgen::collflux::CollinearFluxFactory::get().describe(cf_gr.first).data()));
        coll.emplace_back(&gr);
      }
    plt->draw(coll, "comp_collfluxes", oss.str(), dm);
  }

  return 0;
}
