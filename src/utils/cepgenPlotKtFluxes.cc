/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
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
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/Beam.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/String.h"

int main(int argc, char* argv[]) {
  std::string formfac_type;
  int strfun_type, num_points;
  double kt2, mx;
  bool logy, draw_grid;
  std::string output_file, plotter;

  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("ff,f", "form factors modelling", &formfac_type, "StandardDipole")
      .addOptionalArgument("sf,s", "structure functions modelling", &strfun_type, 301)
      .addOptionalArgument("kt2,k", "parton transverse virtuality (GeV^2)", &kt2, 100.)
      .addOptionalArgument("mx,m", "diffractive state mass (GeV)", &mx, 1.5)
      .addOptionalArgument("npoints,n", "number of x-points to scan", &num_points, 100)
      .addOptionalArgument("output,o", "output file name", &output_file, "flux.scan.output.txt")
      .addOptionalArgument("logy,l", "logarithmic y-scale", &logy, false)
      .addOptionalArgument("draw-grid,g", "draw the x/y grid", &draw_grid, false)
      .addOptionalArgument("plotter,p", "type of plotter to user", &plotter, "")
      .parse();

  cepgen::initialise();
  const double mi = cepgen::PDG::get().mass(cepgen::PDG::proton);
  const double mi2 = mi * mi, mx2 = mx * mx;

  auto ff = cepgen::formfac::FormFactorsFactory::get().build(formfac_type);
  auto sf = cepgen::strfun::StructureFunctionsFactory::get().build(strfun_type);
  std::ofstream out(output_file);
  out << "# form factors: " << ff.get() << "\n"
      << "# structure functions: " << sf.get() << "\n"
      << "# kt2 = " << kt2 << " GeV^2\n"
      << "# mX = " << mx << " GeV\n";
  cepgen::utils::Graph1D graph_el("", "Elastic photon"), graph_inel("", "Inelastic photon"),
      graph_inel_bud("", "Inelastic photon (Budnev)"), graph_glu("", "Gluon (KMR)");
  for (int i = 0; i < num_points; ++i) {
    const double x = i * 1. / num_points;
    const auto f_el = cepgen::Beam::ktFluxNucl(cepgen::Beam::KTFlux::P_Photon_Elastic, x, kt2, ff.get());
    const auto f_inel =
        cepgen::Beam::ktFluxNucl(cepgen::Beam::KTFlux::P_Photon_Inelastic, x, kt2, nullptr, sf.get(), mi2, mx2);
    const auto f_inel_bud =
        cepgen::Beam::ktFluxNucl(cepgen::Beam::KTFlux::P_Photon_Inelastic_Budnev, x, kt2, nullptr, sf.get(), mi2, mx2);
    //const auto f_glu = cepgen::Beam::ktFluxNucl(cepgen::Beam::KTFlux::P_Gluon_KMR, x, kt2, nullptr, nullptr, mi2, mx2);
    out << x << "\t" << f_el << "\t" << f_inel << "\t" << f_inel_bud << /*"\t" << f_glu <<*/ "\n";
    graph_el.addPoint(x, f_el);
    graph_inel.addPoint(x, f_inel);
    graph_inel_bud.addPoint(x, f_inel_bud);
    //graph_glu.addPoint(x, f_glu);
  }
  out.close();
  CG_LOG << "Scan written in \"" << output_file << "\".";

  if (!plotter.empty()) {
    auto plt = cepgen::utils::DrawerFactory::get().build(plotter);
    cepgen::utils::Drawer::Mode dm;
    if (logy)
      dm |= cepgen::utils::Drawer::Mode::logy;
    if (draw_grid)
      dm |= cepgen::utils::Drawer::Mode::grid;

    const auto top_label = cepgen::utils::format("k_{T}^{2} = %g GeV^{2}", kt2) + ", " +
                           cepgen::formfac::FormFactorsFactory::get().describe(formfac_type) + "/" +
                           cepgen::strfun::StructureFunctionsFactory::get().describe(strfun_type);

    graph_el.xAxis().setLabel("#xi");
    graph_el.yAxis().setLabel("#varphi_{T}(#xi, k_{T}^{2})");
    plt->draw({&graph_el, &graph_inel, &graph_inel_bud}, "comp_ktflux", top_label, dm);
  }

  return 0;
}
