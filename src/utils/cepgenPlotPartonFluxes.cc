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
  string formfac_type;
  int strfun_type, num_points;
  double kt2, mx;
  bool logx, logy, draw_grid;
  string output_file, plotter;

  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument(
          "fluxes", "parton fluxe modellings", &fluxes_names, cepgen::PartonFluxFactory::get().modules())
      .addOptionalArgument("ff,f", "form factors modelling", &formfac_type, "StandardDipole")
      .addOptionalArgument("sf,s", "structure functions modelling", &strfun_type, 301)
      .addOptionalArgument("kt2,k", "parton transverse virtuality (GeV^2)", &kt2, 100.)
      .addOptionalArgument("mx,m", "diffractive state mass (GeV)", &mx, 1.5)
      .addOptionalArgument("npoints,n", "number of x-points to scan", &num_points, 100)
      .addOptionalArgument("output,o", "output file name", &output_file, "flux.scan.output.txt")
      .addOptionalArgument("logx", "logarithmic x-scale", &logx, false)
      .addOptionalArgument("logy,l", "logarithmic y-scale", &logy, false)
      .addOptionalArgument("draw-grid,g", "draw the x/y grid", &draw_grid, false)
      .addOptionalArgument("plotter,p", "type of plotter to user", &plotter, "")
      .parse();

  cepgen::initialise();
  const double mx2 = mx * mx;

  ofstream out(output_file);
  out << "# parton fluxes: " << cepgen::utils::repr(fluxes_names) << "\n"
      << "# form factors: " << formfac_type << "\n"
      << "# structure functions: " << strfun_type << "\n"
      << "# kt2 = " << kt2 << " GeV^2\n"
      << "# mX = " << mx << " GeV";
  vector<cepgen::utils::Graph1D> graph_flux;
  vector<std::unique_ptr<cepgen::PartonFlux> > fluxes;
  for (const auto& flux : fluxes_names) {
    auto flux_obj = cepgen::PartonFluxFactory::get().build(
        flux,
        cepgen::ParametersList()
            .set<cepgen::ParametersList>(
                "structureFunctions",
                cepgen::strfun::StructureFunctionsFactory::get().describeParameters(strfun_type).parameters())
            .set<cepgen::ParametersList>(
                "formFactors",
                cepgen::formfac::FormFactorsFactory::get().describeParameters(formfac_type).parameters()));
    graph_flux.emplace_back(flux, flux_obj->name());
    fluxes.emplace_back(std::move(flux_obj));
  }
  for (const auto& x : cepgen::Limits(0., 1. - 1.e-15).generate(num_points)) {
    out << "\n" << x;
    for (size_t i = 0; i < fluxes.size(); ++i) {
      const auto f = (*fluxes.at(i))(x, kt2, mx2);
      out << "\t" << f;
      graph_flux.at(i).addPoint(x, f);
    }
  }
  out.close();
  CG_LOG << "Scan written in \"" << output_file << "\".";

  if (!plotter.empty()) {
    auto plt = cepgen::utils::DrawerFactory::get().build(plotter);
    cepgen::utils::Drawer::Mode dm;
    if (logx)
      dm |= cepgen::utils::Drawer::Mode::logx;
    if (logy)
      dm |= cepgen::utils::Drawer::Mode::logy;
    if (draw_grid)
      dm |= cepgen::utils::Drawer::Mode::grid;

    const auto top_label = cepgen::utils::format("$k_{T}^{2}$ = %g GeV$^{2}$", kt2) + ", " +
                           cepgen::formfac::FormFactorsFactory::get().describe(formfac_type) + "/" +
                           cepgen::strfun::StructureFunctionsFactory::get().describe(strfun_type);

    cepgen::utils::DrawableColl plots;
    for (auto& plt : graph_flux) {
      plt.xAxis().setLabel("$\\xi$");
      plt.yAxis().setLabel("$\\varphi_{T}(\\xi, k_{T}^{2})$");
      plots.emplace_back(&plt);
    }
    plt->draw(plots, "comp_partonflux", top_label, dm);
  }

  return 0;
}
