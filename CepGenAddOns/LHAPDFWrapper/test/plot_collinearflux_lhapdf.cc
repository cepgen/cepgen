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

#include <LHAPDF/LHAPDF.h>

#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/Beam.h"
#include "CepGen/Physics/CollinearFlux.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"

using namespace std;

int main(int argc, char* argv[]) {
  double q2, xmin, xmax;
  string ffmode, set, output, plotter;
  int strfun_type, member, num_points;
  bool logx, logy, draw_grid;
  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("q2", "Virtuality", &q2, 100.)
      .addOptionalArgument("xmin,x", "minimal fractional loss", &xmin, 1.e-5)
      .addOptionalArgument("xmax,X", "maximal fractional loss", &xmax, 1.)
      .addOptionalArgument("formfac,f", "form factors modelling", &ffmode, "StandardDipole")
      .addOptionalArgument("sf,s", "structure functions modelling", &strfun_type, 301)
      .addOptionalArgument("set,s", "PDFset to use", &set, "LUXqed17_plus_PDF4LHC15_nnlo_100")
      .addOptionalArgument("output,o", "Output filename", &output, argv[0])
      .addOptionalArgument("member,m", "PDF member", &member, 0)
      .addOptionalArgument("num-points,n", "Number of points to probe", &num_points, 100)
      .addOptionalArgument("plotter,p", "type of plotter to user", &plotter, "")
      .addOptionalArgument("logx", "logarithmic x-axis", &logx, false)
      .addOptionalArgument("logy,l", "logarithmic y-axis", &logy, false)
      .addOptionalArgument("draw-grid,g", "draw the x/y grid", &draw_grid, false)
      .parse();

  const double lxmin = log10(xmin), lxmax = log10(xmax);

  cepgen::initialise();

  unique_ptr<LHAPDF::PDF> pdf(LHAPDF::mkPDF(set, member));

  //auto sf = cepgen::strfun::StructureFunctionsFactory::get().build(strfun_type);
  auto sf = cepgen::strfun::StructureFunctionsFactory::get().build(
      401, cepgen::ParametersList().set<std::string>("pdfSet", set).set<int>("pdfMember", member));
  auto ff = cepgen::formfac::FormFactorsFactory::get().build(ffmode);

  const cepgen::Limits kt2_limits(0., 1000.);
  const cepgen::CollinearFlux flux(ff.get(), sf.get(), kt2_limits);

  cepgen::utils::Graph1D g_ref, g_cg, g_ratio;
  for (int i = 0; i < num_points; ++i) {
    const double x =
        (!logx) ? xmin + i * (xmax - xmin) / (num_points - 1) : pow(10, lxmin + i * (lxmax - lxmin) / (num_points - 1));
    const double xfx = pdf->xfxQ2(22, x, q2);
    const double pdf = flux(x, 0.938, cepgen::Beam::KTFlux::P_Photon_Elastic_Budnev);
    cout << x << "\t" << xfx << "\t" << pdf << "\t" << pdf / xfx << endl;
    g_ref.addPoint(x, xfx);
    g_cg.addPoint(x, pdf);
    g_ratio.addPoint(x, pdf / xfx);
  }

  if (!plotter.empty()) {
    auto plt = cepgen::utils::DrawerFactory::get().build(plotter);
    cepgen::utils::Drawer::Mode dm;
    if (logx)
      dm |= cepgen::utils::Drawer::Mode::logx;
    if (logy)
      dm |= cepgen::utils::Drawer::Mode::logy;
    if (draw_grid)
      dm |= cepgen::utils::Drawer::Mode::grid;
    cepgen::utils::DrawableColl mg;
    for (auto* gr : {&g_ref, &g_cg, &g_ratio}) {
      gr->xAxis().setLabel("$x$");
      gr->yAxis().setLabel("$f_{\\gamma}(x)$");
      mg.emplace_back(gr);
    }
    plt->draw(mg, output);
  }
  return 0;
}
