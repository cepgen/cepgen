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
#include "CepGen/Generator.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/Utils.h"
#include "CepGen/StructureFunctions/SigmaRatio.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/String.h"

using namespace std;

int main(int argc, char* argv[]) {
  vector<int> sigrat_types;
  double q2, w2;
  cepgen::Limits x_range;
  int var, num_points;
  string output_file, plotter;
  bool logx, logy, draw_grid;

  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("sr,s", "longitud./transv. cross section ratio modelling", &sigrat_types)
      .addOptionalArgument("q2,q", "parton virtuality (GeV^2)", &q2, -1.)
      .addOptionalArgument("w2,w", "scattered particle squared mass (GeV^2/c^4)", &w2, -1.)
      .addOptionalArgument("var,t", "variable to study (0=xBj, 1=w)", &var, 0)
      .addOptionalArgument("xrange,x", "Bjorken x range", &x_range, cepgen::Limits{1.e-7, 1.})
      .addOptionalArgument("npoints,n", "number of x-points to scan", &num_points, 500)
      .addOptionalArgument("output,o", "output file name", &output_file, "strfuns.scan.output.txt")
      .addOptionalArgument("logx", "logarithmic x-axis", &logx, false)
      .addOptionalArgument("logy,l", "logarithmic y-axis", &logy, false)
      .addOptionalArgument("draw-grid,g", "draw the x/y grid", &draw_grid, false)
      .addOptionalArgument("plotter,p", "type of plotter to user", &plotter, "")
      .parse();

  if (q2 < 0. && w2 < 0.)
    throw CG_FATAL("main") << "Either a Q^2 or a w^2 must be provided!";

  cepgen::initialise();

  string var_name, var_unit;
  switch (var) {
    case 0:
      var_name = "x_{Bj}";
      break;
    case 1:
      var_name = "w";
      var_unit = "GeV";
      break;
    case 2:
      var_name = "w^{2}";
      var_unit = "GeV$^{2}$";
      break;
    default:
      throw CG_FATAL("main") << "Unsupported variable to be plotted!";
  }

  string fixed_var;
  if (q2 > 0.)
    fixed_var = "$Q^{2}$";
  else if (w2 > 0.)
    fixed_var = "$w^{2}$";

  ofstream out(output_file);
  out << "# sigma ratios: ";
  string sep;
  for (const auto& sr_type : sigrat_types)
    out << sep << sr_type, sep = ", ";
  out << "\n"
      << "# x in [" << x_range << "]\n";

  const float mp = cepgen::PDG::get().mass(2212), mp2 = mp * mp;

  vector<unique_ptr<cepgen::sigrat::Parameterisation> > sigrats;
  vector<cepgen::utils::Graph1D> g_sigrats;
  if (sigrat_types.empty())
    sigrat_types = cepgen::SigmaRatiosFactory::get().modules();
  for (const auto& sr_type : sigrat_types) {
    auto sr = cepgen::SigmaRatiosFactory::get().build(sr_type);
    const auto sr_name = cepgen::SigmaRatiosFactory::get().describe(sr_type);
    g_sigrats.emplace_back(to_string(sr_type), sr_name);
    sigrats.emplace_back(move(sr));
  }
  const auto xvals = x_range.generate(num_points, logx);
  for (int i = 0; i < num_points; ++i) {
    const auto& x = xvals.at(i);
    out << x << "\t";
    size_t j = 0;
    for (auto& sr : sigrats) {
      double xbj;
      switch (var) {
        case 0:
          xbj = x;
          break;
        case 1:
          xbj = cepgen::utils::xBj(q2, mp2, x * x);
          break;
        case 2:
          xbj = cepgen::utils::xBj(q2, mp2, x);
          break;
        default:
          xbj = 0.;
          break;
      }
      double err;
      const auto sigrat = (*sr)(xbj, q2, err);
      out << "\t" << sigrat;
      g_sigrats.at(j).addPoint(x, sigrat);
      ++j;
    }
    out << "\n";
  }
  CG_LOG << "Scan written in \"" << output_file << "\".";
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

    cepgen::utils::DrawableColl mg;
    for (auto& gr : g_sigrats) {
      gr.xAxis().setLabel("$" + var_name + "$" + (!var_unit.empty() ? " (" + var_unit + ")" : ""));
      gr.yAxis().setLabel("$\\sigma_{L}/\\sigma_{T} = R(" + var_name + ", Q^{2})$");
      mg.emplace_back(&gr);
    }
    plt->draw(mg, "comp_sigrat", fixed_var + cepgen::utils::format(" = %g GeV$^{2}$", q2 > 0. ? q2 : w2), dm);
  }

  return 0;
}
