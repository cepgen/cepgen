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

#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Modules/FormFactorsFactory.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/Message.h"
#include "CepGen/Utils/String.h"

using namespace std;

int main(int argc, char* argv[]) {
  int num_points;
  string output_file, plotter;
  bool logx, logy, draw_grid;
  cepgen::Limits q2range, yrange;
  vector<string> modules;

  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("modules,m", "types of form factors", &modules, cepgen::FormFactorsFactory::get().modules())
      .addOptionalArgument("q2range,q", "parton virtuality range (GeV^2)", &q2range, cepgen::Limits{1., 2.5})
      .addOptionalArgument("yrange,y", "y range", &yrange)
      .addOptionalArgument("npoints,n", "number of x-points to scan", &num_points, 500)
      .addOptionalArgument("output,o", "output file name", &output_file, "formfacs.scan.output.txt")
      .addOptionalArgument("logx", "logarithmic x-axis", &logx, false)
      .addOptionalArgument("logy,l", "logarithmic y-axis", &logy, false)
      .addOptionalArgument("draw-grid,g", "draw the x/y grid", &draw_grid, false)
      .addOptionalArgument("plotter,p", "type of plotter to user", &plotter, "")
      .parse();

  cepgen::initialise();

  ofstream out(output_file);
  out << "# form factors: " << cepgen::utils::merge(modules, ",");

  out << "\n"
      << "# q2 range: " << q2range << " GeV^2\n";

  vector<unique_ptr<cepgen::formfac::Parameterisation> > form_factors;
  vector<cepgen::utils::Graph1D> g_form_factors_fe, g_form_factors_fm;
  for (const auto& ff_type : modules) {
    form_factors.emplace_back(cepgen::FormFactorsFactory::get().build(ff_type));
    const auto ff_desc = cepgen::FormFactorsFactory::get().describe(ff_type);
    g_form_factors_fe.emplace_back("fe_" + ff_type, ff_desc);
    g_form_factors_fm.emplace_back("fm_" + ff_type, ff_desc);
  }
  for (const auto& q2 : q2range.generate(num_points, logx)) {
    out << q2 << "\t";
    size_t j = 0;
    for (auto& ff : form_factors) {
      const auto form_factor = (*ff)(q2);
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
    if (logx)
      dm |= cepgen::utils::Drawer::Mode::logx;
    if (logy)
      dm |= cepgen::utils::Drawer::Mode::logy;
    if (draw_grid)
      dm |= cepgen::utils::Drawer::Mode::grid;

    for (auto& canv : map<pair<string, string>, vector<cepgen::utils::Graph1D> >{
             {{"fe", "$F_{E}$"}, g_form_factors_fe}, {{"fm", "$F_{M}$"}, g_form_factors_fm}}) {
      cepgen::utils::DrawableColl mp;
      for (auto& gr : canv.second) {
        gr.xAxis().setLabel("Q$^{2}$ (GeV$^{2}$)");
        gr.yAxis().setLabel(canv.first.second);
        if (yrange.valid())
          gr.yAxis().setRange(yrange);
        mp.emplace_back(&gr);
      }
      plt->draw(mp, "comp_" + canv.first.first, "", dm);
    }
  }

  return 0;
}
