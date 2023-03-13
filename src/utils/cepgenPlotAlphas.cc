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
#include "CepGen/Modules/CouplingFactory.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Physics/Coupling.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/String.h"

using namespace std;

int main(int argc, char* argv[]) {
  cepgen::Limits q_range;
  int num_points;
  string output_file, plotter;
  bool q2mode, logx, logy, draw_grid;

  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("qrange,q", "virtuality range (GeV)", &q_range, cepgen::Limits{1., 101.})
      .addOptionalArgument("q2mode", "plot as a function of Q^2", &q2mode, false)
      .addOptionalArgument("npoints,n", "number of x-points to scan", &num_points, 100)
      .addOptionalArgument("output,o", "output file name", &output_file, "alphas.scan.output.txt")
      .addOptionalArgument("logx", "logarithmic x-scale", &logx, false)
      .addOptionalArgument("logy,l", "logarithmic y-scale", &logy, false)
      .addOptionalArgument("draw-grid,g", "draw the x/y grid", &draw_grid, false)
      .addOptionalArgument("plotter,p", "type of plotter to user", &plotter, "")
      .parse();

  cepgen::initialise();

  struct alpha_t {
    string name;
    vector<double> vals;
    cepgen::utils::Graph1D graph;
  };
  vector<alpha_t> alphas, alphaem;

  const auto qvals = q_range.generate(num_points, logx);

  {
    // alphaS(Q) modellings part
    size_t i = 0;
    for (const auto& mod : cepgen::AlphaSFactory::get().modules()) {
      const auto& algo = cepgen::AlphaSFactory::get().build(
          mod /*, cepgen::ParametersList().set<double>("asmur", 0.35).set<double>("mur", 1.4142)*/);
      alphas.emplace_back(alpha_t{
          mod,
          vector<double>(num_points),
          cepgen::utils::Graph1D(
              mod, cepgen::utils::replace_all(cepgen::AlphaSFactory::get().describe(mod), "alpha(S)", "\\alpha_{S}"))});
      auto& as = alphas[i++];
      for (size_t j = 0; j < qvals.size(); ++j) {
        const auto val = (*algo)(qvals[j]);
        as.vals[j] = val;
        as.graph.addPoint(q2mode ? qvals[j] * qvals[j] : qvals[j], val);
      }
    }
    // alphaEM(Q) modellings part
    i = 0;
    for (const auto& mod : cepgen::AlphaEMFactory::get().modules()) {
      const auto& algo = cepgen::AlphaEMFactory::get().build(mod);
      alphaem.emplace_back(alpha_t{
          mod,
          vector<double>(num_points),
          cepgen::utils::Graph1D(
              mod,
              cepgen::utils::replace_all(cepgen::AlphaEMFactory::get().describe(mod), "alpha(EM)", "\\alpha_{EM}"))});
      auto& aem = alphaem[i++];
      for (size_t j = 0; j < qvals.size(); ++j) {
        const auto val = (*algo)(qvals[j]);
        aem.vals[j] = val;
        aem.graph.addPoint(q2mode ? qvals[j] * qvals[j] : qvals[j], val);
      }
    }
  }

  // output ascii file
  ofstream out(output_file);
  out << "#";
  for (const auto& smp : alphas)
    out << "\t" << smp.name;
  for (const auto& smp : alphaem)
    out << "\t" << smp.name;
  for (size_t i = 0; i < qvals.size(); ++i) {
    out << "\n" << (q2mode ? qvals[i] * qvals[i] : qvals[i]);
    for (const auto& smp : alphas)
      out << "\t" << smp.vals[i];
    for (const auto& smp : alphaem)
      out << "\t" << smp.vals[i];
  }

  // drawing part

  if (!plotter.empty()) {
    auto plt = cepgen::DrawerFactory::get().build(plotter);
    cepgen::utils::Drawer::Mode dm;
    if (logx)
      dm |= cepgen::utils::Drawer::Mode::logx;
    if (logy)
      dm |= cepgen::utils::Drawer::Mode::logy;
    if (draw_grid)
      dm |= cepgen::utils::Drawer::Mode::grid;
    string xlabel = q2mode ? "Q^{2} (GeV^{2})" : "Q (GeV)", spectrum = q2mode ? "Q^{2}" : "Q";

    {
      cepgen::utils::DrawableColl mp;
      for (size_t i = 0; i < alphas.size(); ++i) {
        alphas[i].graph.xAxis().setLabel(xlabel);
        alphas[i].graph.yAxis().setLabel("$\\alpha_{S}(" + spectrum + ")$");
        mp.emplace_back(&alphas[i].graph);
        //const auto descr = cepgen::utils::replace_all(cepgen::AlphaSFactory::get().describe(alphas[i].name),
        //                                              {{" alphaS", ""}, {" evolution algorithm", ""}});
      }
      plt->draw(mp, "comp_alphas", cepgen::utils::s("$\\alpha_{S}$ modelling", alphas.size(), false), dm);
    }
    {
      cepgen::utils::DrawableColl mp;
      for (size_t i = 0; i < alphaem.size(); ++i) {
        alphaem[i].graph.xAxis().setLabel(xlabel);
        alphaem[i].graph.yAxis().setLabel("$\\alpha_{EM}$(" + spectrum + ")");
        mp.emplace_back(&alphaem[i].graph);
        //const auto descr = cepgen::utils::replace_all(cepgen::AlphaEMFactory::get().describe(alphaem[i].name),
        //                                              {{" alphaS", ""}, {" evolution algorithm", ""}});
      }
      plt->draw(mp, "comp_alphaem", cepgen::utils::s("$\\alpha_{EM}$ modelling", alphaem.size(), false), dm);
    }
  }
  return 0;
}
