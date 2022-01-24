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
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/String.h"

using namespace std;

int main(int argc, char* argv[]) {
  vector<int> strfun_types;
  double q2, xmin, xmax, ymin, ymax;
  int var, num_points;
  string output_file, plotter;
  bool logx, logy, draw_grid;

  const double kInvalid = -999.999;

  cepgen::ArgumentsParser(argc, argv)
      .addArgument("sf,s", "structure functions modelling", &strfun_types)
      .addOptionalArgument("q2,q", "parton virtuality (GeV^2)", &q2, 10.)
      .addOptionalArgument("var,t", "variable to study (0=xBj, 1=w)", &var, 0)
      .addOptionalArgument("xmax,m", "minimal Bjorken x", &xmin, 1.e-7)
      .addOptionalArgument("xmax,M", "maximal Bjorken x", &xmax, 1.)
      .addOptionalArgument("ymin", "minimum plotting range for y", &ymin, kInvalid)
      .addOptionalArgument("ymax", "maximum plotting range for y", &ymax, kInvalid)
      .addOptionalArgument("npoints,n", "number of x-points to scan", &num_points, 500)
      .addOptionalArgument("output,o", "output file name", &output_file, "strfuns.scan.output.txt")
      .addOptionalArgument("logx", "logarithmic x-axis", &logx, false)
      .addOptionalArgument("logy,l", "logarithmic y-axis", &logy, false)
      .addOptionalArgument("draw-grid,g", "draw the x/y grid", &draw_grid, false)
      .addOptionalArgument("plotter,p", "type of plotter to user", &plotter, "")
      .parse();

  const double lxmin = log10(xmin), lxmax = log10(xmax);

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
      var_unit = "GeV^{2}";
      break;
    default:
      throw CG_FATAL("main") << "Unsupported variable to be plotted!";
  }

  ofstream out(output_file);
  out << "# structure functions: ";
  string sep;
  for (const auto& sf_type : strfun_types)
    out << sep << sf_type, sep = ", ";
  out << "\n"
      << "# x in [" << xmin << ", " << xmax << "]\n";

  const float mp = cepgen::PDG::get().mass(2212), mp2 = mp * mp;

  cepgen::utils::Drawer::Mode dm;
  if (logx)
    dm |= cepgen::utils::Drawer::Mode::logx;
  if (logy)
    dm |= cepgen::utils::Drawer::Mode::logy;
  if (draw_grid)
    dm |= cepgen::utils::Drawer::Mode::grid;

  vector<unique_ptr<cepgen::strfun::Parameterisation> > strfuns;
  vector<cepgen::utils::Graph1D> g_strfuns_f2, g_strfuns_fl, g_strfuns_fe, g_strfuns_fm, g_strfuns_w1, g_strfuns_w2;
  for (const auto& sf_type : strfun_types) {
    auto sf = cepgen::strfun::StructureFunctionsFactory::get().build(sf_type);
    const auto sf_name = cepgen::strfun::StructureFunctionsFactory::get().describe(sf_type);
    g_strfuns_f2.emplace_back("f2_" + sf_type, sf_name);
    g_strfuns_fl.emplace_back("fl_" + sf_type, sf_name);
    g_strfuns_fe.emplace_back("fe_" + sf_type, sf_name);
    g_strfuns_fm.emplace_back("fm_" + sf_type, sf_name);
    g_strfuns_w1.emplace_back("w1_" + sf_type, sf_name);
    g_strfuns_w2.emplace_back("w2_" + sf_type, sf_name);
    strfuns.emplace_back(move(sf));
  }
  for (int i = 0; i < num_points; ++i) {
    const double x =
        (!logx) ? xmin + i * (xmax - xmin) / (num_points - 1) : pow(10, lxmin + i * (lxmax - lxmin) / (num_points - 1));
    out << x << "\t";
    size_t j = 0;
    for (auto& sf : strfuns) {
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
      out << "\t" << sf->F2(xbj, q2) << "\t" << sf->FL(xbj, q2);
      g_strfuns_f2.at(j).addPoint(x, sf->F2(xbj, q2));
      g_strfuns_fl.at(j).addPoint(x, sf->FL(xbj, q2));
      g_strfuns_fe.at(j).addPoint(x, sf->FE(xbj, q2));
      g_strfuns_fm.at(j).addPoint(x, sf->FM(xbj, q2));
      g_strfuns_w1.at(j).addPoint(x, sf->W1(xbj, q2));
      g_strfuns_w2.at(j).addPoint(x, sf->W2(xbj, q2));
      ++j;
    }
    out << "\n";
  }
  CG_LOG << "Scan written in \"" << output_file << "\".";
  out.close();

  if (!plotter.empty()) {
    auto plt = cepgen::utils::DrawerFactory::get().build(plotter);
    for (auto& canv : map<pair<string, string>, vector<cepgen::utils::Graph1D> >{{{"f2", "F_{2}"}, g_strfuns_f2},
                                                                                 {{"fl", "F_{L}"}, g_strfuns_fl},
                                                                                 {{"fe", "F_{E}"}, g_strfuns_fe},
                                                                                 {{"fm", "F_{M}"}, g_strfuns_fm},
                                                                                 {{"w1", "W_{1}"}, g_strfuns_w1},
                                                                                 {{"w2", "W_{2}"}, g_strfuns_w2}}) {
      cepgen::utils::DrawableColl mp;
      for (auto& p : canv.second) {
        p.xAxis().setLabel(var_name + (!var_unit.empty() ? " (" + var_unit + ")" : ""));
        p.yAxis().setLabel(canv.first.second + (!var_name.empty() ? "(" + var_name + ", Q^{2})" : ""));
        if (ymin != kInvalid)
          p.yAxis().setMinimum(ymin);
        if (ymax != kInvalid)
          p.yAxis().setMaximum(ymax);
        mp.emplace_back(&p);
      }
      plt->draw(mp, "sfcomp_" + canv.first.first, cepgen::utils::format("Q^{2} = %g GeV^{2}", q2), dm);
    }
  }

  return 0;
}
