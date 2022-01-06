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

#include "CepGen/Cards/Handler.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ExportModule.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/CardsHandlerFactory.h"
#include "CepGen/Parameters.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Processes/Process.h"
#include "CepGen/Utils/AbortHandler.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Logger.h"
#include "CepGen/Utils/String.h"

using namespace std;

int main(int argc, char* argv[]) {
  string input_config, output_file, scan;
  int npoints;
  double min_value, max_value;
  vector<double> points;

  cepgen::ArgumentsParser parser(argc, argv);
  parser.addArgument("config,i", "base configuration", &input_config)
      .addOptionalArgument("scan,s", "type of scan to perform", &scan, "ptmin")
      .addOptionalArgument("min,l", "minimum value of scan", &min_value, 1.)
      .addOptionalArgument("max,H", "maximum value of scan", &max_value, 11.)
      .addOptionalArgument("num-points,n", "number of points to consider", &npoints, 10)
      .addOptionalArgument("points,p", "list of points to consider", &points, vector<double>{})
      .addOptionalArgument("output,o", "output file", &output_file, "xsect.dat")
      .parse();

  cepgen::Generator mg;
  mg.setParameters(cepgen::card::Handler::parse(input_config));

  if (!parser.extra_config().empty())
    mg.setParameters(cepgen::card::CardsHandlerFactory::get()
                         .build(cepgen::card::gCommandLineHandler,
                                cepgen::ParametersList().set<std::vector<std::string> >("args", parser.extra_config()))
                         ->parse("", mg.parametersPtr()));

  CG_LOG << mg.parameters();

  ofstream xsect_file(output_file);
  if (!xsect_file.is_open())
    throw CG_FATAL("main") << "Output file \"" << output_file << "\" cannot be opened!";
  xsect_file << "# " << scan << "\txsect (pb)\td(xsect) (pb)\n";

  auto& par = mg.parametersRef();
  //--- ensure nothing is written in the output sequence
  par.outputModulesSequence().clear();

  if (points.empty())
    for (int i = 0; i <= npoints; ++i)
      points.emplace_back(min_value + (max_value - min_value) * i / npoints);

  cepgen::utils::AbortHandler();

  auto& kin = par.process().kinematics();
  double cross_section, err_cross_section;
  for (const auto& value : points) {
    try {
      if (scan == "sqrtS") {
        kin.incomingBeams().setSqrtS(value);
      } else if (scan == "abseta") {
        kin.cuts().central.eta_single().min() = -value;
        kin.cuts().central.eta_single().max() = +value;
      } else if (scan == "absrap") {
        kin.cuts().central.rapidity_single().min() = -value;
        kin.cuts().central.rapidity_single().max() = +value;
      } else if (scan == "mpart") {
        auto prop = cepgen::PDG::get()(par.process().event()[cepgen::Particle::CentralSystem][0].pdgId());
        prop.mass = value;
        cepgen::PDG::get().define(prop);
        par.process().clear();
      } else {
        auto modif = cepgen::ParametersList().set<double>(scan, value);
        kin.setParameters(modif);
        CG_LOG << modif << "\n\n" << kin.cuts();
      }
      CG_LOG << "Scan of \"" << scan << "\". Value = " << value << ".";
      mg.computeXsection(cross_section, err_cross_section);
      string out_line = cepgen::utils::format("%.2f\t%.8e\t%.8e\n", value, cross_section, err_cross_section);
      xsect_file << out_line;
      CG_LOG << out_line;
      xsect_file.flush();
    } catch (const cepgen::utils::RunAbortedException&) {
      CG_LOG << "Run aborted!";
      break;
    }
  }

  return 0;
}
