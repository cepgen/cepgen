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

#include "CepGen/Core/Exception.h"
#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/String.h"

using namespace std;

int main(int argc, char* argv[]) {
  int mode, strfun_type, num_points;
  double mx, q2min, q2max;
  string output_file, plotter;
  bool logy, draw_grid;

  cepgen::ArgumentsParser(argc, argv)
      .addArgument("mode,t", "beam modelling", &mode, (int)cepgen::Beam::Mode::ProtonElastic)
      .addOptionalArgument("mx,x", "diffractive mass (GeV/c^2)", &mx, 100.)
      .addOptionalArgument("sf,s", "structure functions modelling", &strfun_type, 11)
      .addOptionalArgument("q2min,m", "minimal parton virtuality (GeV^2)", &q2min, 1.)
      .addOptionalArgument("q2max,M", "maximal parton virtuality (GeV^2)", &q2max, 2.5)
      .addOptionalArgument("npoints,n", "number of x-points to scan", &num_points, 500)
      .addOptionalArgument("output,o", "output file name", &output_file, "formfacs.scan.output.txt")
      .addOptionalArgument("logy,l", "logarithmic y-axis", &logy, false)
      .addOptionalArgument("draw-grid,g", "draw the x/y grid", &draw_grid, false)
      .addOptionalArgument("plotter,p", "type of plotter to user", &plotter, "")
      .parse();

  cepgen::initialise();

  ofstream out(output_file);
  out << "# form factors: ";
  string sep;
  for (const auto& ff_type : cepgen::formfac::FormFactorsFactory::get().modules())
    out << sep << ff_type, sep = ", ";

  auto sf = cepgen::strfun::StructureFunctionsFactory::get().build(strfun_type);
  out << "\n"
      << "# structure functions: " << sf.get() << "\n"
      << "# q2 in [" << q2min << ", " << q2max << "] GeV^2\n";

  vector<unique_ptr<cepgen::formfac::Parameterisation> > form_factors;
  vector<cepgen::utils::Graph1D> g_form_factors_fe, g_form_factors_fm;
  for (const auto& ff_type : cepgen::formfac::FormFactorsFactory::get().modules()) {
    form_factors.emplace_back(cepgen::formfac::FormFactorsFactory::get().build(ff_type));
    const auto ff_desc = cepgen::formfac::FormFactorsFactory::get().describe(ff_type);
    g_form_factors_fe.emplace_back("fe_" + ff_type, ff_desc);
    g_form_factors_fm.emplace_back("fm_" + ff_type, ff_desc);
  }
  for (int i = 0; i < num_points; ++i) {
    const double q2 = q2min + i * (q2max - q2min) / (num_points - 1);
    out << q2 << "\t";
    size_t j = 0;
    for (auto& ff : form_factors) {
      const auto form_factor = (*ff)((cepgen::Beam::Mode)mode, q2, mx, sf.get());
      out << "\t" << form_factor.FE << "\t" << form_factor.FM;
      g_form_factors_fe.at(j).addPoint(q2, form_factor.FE);
      g_form_factors_fm.at(j).addPoint(q2, form_factor.FM);
      ++j;
    }
    out << "\n";
  }
  CG_LOG << "Scan written in \"" << output_file << "\".";
  out.close();

  if (!plotter.empty()) {
    auto plt = cepgen::utils::DrawerFactory::get().build(plotter);
    cepgen::utils::Drawer::Mode dm;
    if (logy)
      dm |= cepgen::utils::Drawer::Mode::logy;
    if (draw_grid)
      dm |= cepgen::utils::Drawer::Mode::grid;

    for (auto& canv : map<pair<string, string>, vector<cepgen::utils::Graph1D> >{
             {{"fe", "F_{E}"}, g_form_factors_fe}, {{"fm", "F_{M}"}, g_form_factors_fm}}) {
      cepgen::utils::DrawableColl mp;
      for (auto& gr : canv.second) {
        gr.xAxis().setLabel("Q^{2} (GeV^{2})");
        gr.yAxis().setLabel(canv.first.second);
        mp.emplace_back(&gr);
      }
      plt->draw(mp, "comp_" + canv.first.first, cepgen::utils::format("M_{X} = %g GeV/c^{2}", mx), dm);
    }
  }

  return 0;
}
