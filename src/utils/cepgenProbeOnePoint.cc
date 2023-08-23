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

#include "CepGen/Cards/Handler.h"
#include "CepGen/Generator.h"
#include "CepGen/Parameters.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Message.h"

using namespace std;

int main(int argc, char* argv[]) {
  string input_card;
  vector<double> point;
  bool enable_plugins, verbose;

  cepgen::ArgumentsParser(argc, argv)
      .addArgument("input,i", "input card", &input_card)
      .addOptionalArgument("point,p", "point to test", &point, vector<double>(12, 0.3))
      .addOptionalArgument("verbose,v", "high verbosity mode", &verbose, false)
      .addOptionalArgument("enable-plugins,m", "enable the external plugins", &enable_plugins, false)
      .parse();

  cepgen::Generator gen;
  gen.setParameters(cepgen::card::Handler::parseFile(input_card));

  const auto ndim = gen.parameters()->process().ndim();
  if (point.size() < 2)
    point = vector<double>(ndim, point[0]);
  else if (point.size() != ndim)
    point.resize(ndim);

  if (verbose)
    CG_LOG_LEVEL(debugInsideLoop);

  if (!enable_plugins) {
    gen.parametersPtr()->clearEventModifiersSequence();
    gen.parametersPtr()->clearEventExportersSequence();
  }

  CG_LOG << gen.parameters() << "\n\t"
         << "point: " << point;
  const double weight = gen.computePoint(point);
  CG_LOG << "weight: " << weight;

  return 0;
}
