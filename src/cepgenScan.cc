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

#include "CepGen/Cards/Handler.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/RunParameters.h"
#include "CepGen/EventFilter/EventExporter.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/CardsHandlerFactory.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/AbortHandler.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Graph.h"
#include "CepGen/Utils/Logger.h"
#include "CepGen/Utils/String.h"

using namespace std;

int main(int argc, char* argv[]) {
  string input_config, output_file, scan, plotter;
  int npoints;
  double min_value, max_value;
  vector<double> points;
  bool draw_grid, logy;

  cepgen::ArgumentsParser parser(argc, argv);
  parser.addArgument("config,i", "base configuration", &input_config)
      .addOptionalArgument("scan,s", "type of scan to perform", &scan, "ptmin")
      .addOptionalArgument("min,m", "minimum value of scan", &min_value, 1.)
      .addOptionalArgument("max,M", "maximum value of scan", &max_value, 11.)
      .addOptionalArgument("num-points,n", "number of points to consider", &npoints, 10)
      .addOptionalArgument("points", "list of points to consider", &points, vector<double>{})
      .addOptionalArgument("output,o", "output file", &output_file, "xsect.dat")
      .addOptionalArgument("logy,l", "logarithmic y-scale", &logy, false)
      .addOptionalArgument("draw-grid,g", "draw the x/y grid", &draw_grid, false)
      .addOptionalArgument("plotter,p", "type of plotter to user", &plotter, "")
      .parse();

  cepgen::Generator mg;
  mg.setRunParameters(cepgen::card::Handler::parseFile(input_config));

  if (!parser.extra_config().empty())
    mg.setRunParameters(cepgen::CardsHandlerFactory::get()
                            .build(".cmd", cepgen::ParametersList().set<vector<string> >("args", parser.extra_config()))
                            ->parseString("", &mg.runParameters()));

  CG_LOG << mg.runParameters();

  ofstream xsect_file(output_file);
  if (!xsect_file.is_open())
    throw CG_FATAL("main") << "Output file \"" << output_file << "\" cannot be opened!";
  xsect_file << "# " << scan << "\txsect (pb)\td(xsect) (pb)\n";

  auto& par = mg.runParameters();
  //--- ensure nothing is written in the output sequence
  par.eventExportersSequence().clear();

  if (points.empty())
    for (int i = 0; i <= npoints; ++i)
      points.emplace_back(min_value + (max_value - min_value) * i / npoints);

  cepgen::utils::AbortHandler();

  cepgen::utils::Graph1D graph("comp_sigma_gen");

  auto& kin = par.process().kinematics();
  string scan_str = scan;
  for (const auto& value : points) {
    try {
      if (scan == "sqrtS") {
        kin.incomingBeams().setSqrtS(value);
        scan_str = "$\\sqrt{s}$ (GeV)";
      } else if (scan == "abseta") {
        kin.cuts().central.eta_single.min() = -value;
        kin.cuts().central.eta_single.max() = +value;
        scan_str = "$|\\eta|$";
      } else if (scan == "absrap") {
        kin.cuts().central.rapidity_single.min() = -value;
        kin.cuts().central.rapidity_single.max() = +value;
        scan_str = "$|y|$";
      } else if (cepgen::utils::startsWith(scan, "m:")) {
        const auto tok = cepgen::utils::split(scan, ':');
        if (tok.size() > 2)
          throw CG_FATAL("main") << "Invalid mass scan defined: should follow the \"m:<pdgid int>\" convention!";
        const cepgen::pdgid_t pdg = stod(tok.at(1));
        cepgen::PDG::get()[pdg].mass = value;
        scan_str = "$m_{" + cepgen::PDG::get()(pdg).name + "}$ (GeV)";
      } else {
        auto modif = cepgen::ParametersList().set<double>(scan, value);
        kin.setParameters(modif);
      }
      CG_LOG << "Scan of \"" << scan << "\". Value = " << value << ".";
      const auto cross_section = mg.computeXsection();
      string out_line = cepgen::utils::format("%.2f\t%.8e\t%.8e\n", value, cross_section, cross_section.uncertainty());
      graph.addPoint(value, cross_section, 0., cross_section.uncertainty());
      xsect_file << out_line;
      CG_LOG << out_line;
      xsect_file.flush();
    } catch (const cepgen::utils::RunAbortedException&) {
      CG_LOG << "Run aborted!";
      break;
    }
  }

  if (!plotter.empty()) {
    cepgen::utils::Drawer::Mode dm;
    if (logy)
      dm |= cepgen::utils::Drawer::Mode::logy;
    if (draw_grid)
      dm |= cepgen::utils::Drawer::Mode::grid;
    graph.xAxis().setLabel(scan_str);
    graph.yAxis().setLabel("$\\sigma_{gen}$ (pb)");
    auto plt = cepgen::DrawerFactory::get().build(plotter);
    plt->draw(graph, dm);
  }

  return 0;
}
