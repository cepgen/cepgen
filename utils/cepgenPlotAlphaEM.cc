/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2025  Laurent Forthomme
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
  cepgen::initialise();

  cepgen::Limits q_range;
  int num_points;
  string output_file, plotter;
  vector<string> models;
  bool q2mode, logx, logy, draw_grid, ratio_plot;

  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("models,m", "models to draw", &models, cepgen::AlphaEMFactory::get().modules())
      .addOptionalArgument("qrange,q", "virtuality range (GeV)", &q_range, cepgen::Limits{1., 101.})
      .addOptionalArgument("q2mode", "plot as a function of Q^2", &q2mode, false)
      .addOptionalArgument("npoints,n", "number of x-points to scan", &num_points, 100)
      .addOptionalArgument("output,o", "output file name", &output_file, "alphaem.scan.output.txt")
      .addOptionalArgument("logx", "logarithmic x-scale", &logx, false)
      .addOptionalArgument("logy,l", "logarithmic y-scale", &logy, false)
      .addOptionalArgument("draw-grid,g", "draw the x/y grid", &draw_grid, false)
      .addOptionalArgument("ratio,r", "draw the ratio plot", &ratio_plot, false)
      .addOptionalArgument("plotter,p", "type of plotter to user", &plotter, "")
      .parse();

  struct alpha_t {
    string name;
    vector<double> vals;
    cepgen::utils::Graph1D graph;
  };
  vector<alpha_t> alphaem;

  const auto q_values = q_range.generate(num_points, logx);
  {  // alphaEM(Q) modellings part
    size_t i = 0;
    for (const auto& mod : models) {
      const auto algo = cepgen::AlphaEMFactory::get().build(mod);
      alphaem.emplace_back(alpha_t{
          mod,
          vector<double>(num_points),
          cepgen::utils::Graph1D(
              mod,
              cepgen::utils::replaceAll(cepgen::AlphaEMFactory::get().describe(mod), "alpha(EM)", "\\alpha_{EM}"))});
      auto& aem = alphaem[i++];
      for (size_t j = 0; j < q_values.size(); ++j) {
        const auto val = (*algo)(q_values[j]);
        aem.vals[j] = val;
        aem.graph.addPoint(q2mode ? q_values[j] * q_values[j] : q_values[j], val);
      }
    }
  }

  // output ascii file
  ofstream out(output_file);
  out << "#";
  for (const auto& smp : alphaem)
    out << "\t" << smp.name;
  for (size_t i = 0; i < q_values.size(); ++i) {
    out << "\n" << (q2mode ? q_values[i] * q_values[i] : q_values[i]);
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
    if (ratio_plot)
      dm |= cepgen::utils::Drawer::Mode::ratio;
    string xlabel = q2mode ? "Q^{2} (GeV^{2})" : "Q (GeV)", spectrum = q2mode ? "Q^{2}" : "Q";

    {
      cepgen::utils::DrawableColl mp;
      for (size_t i = 0; i < alphaem.size(); ++i) {
        alphaem[i].graph.xAxis().setLabel(xlabel);
        alphaem[i].graph.yAxis().setLabel("$\\alpha_{EM}$(" + spectrum + ")");
        mp.emplace_back(&alphaem[i].graph);
      }
      plt->draw(mp, "comp_alphaem", cepgen::utils::s("$\\alpha_{EM}$ modelling", alphaem.size(), false), dm);
    }
  }
  return 0;
}
