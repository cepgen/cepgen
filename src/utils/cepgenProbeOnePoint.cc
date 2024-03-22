/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2024  Laurent Forthomme
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

#include "CepGen/Core/RunParameters.h"
#include "CepGen/Generator.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Message.h"

using namespace std;

int main(int argc, char* argv[]) {
  string input_card;
  int num_scans;
  using point_t = vector<double>;
  point_t user_point;

  const auto args = cepgen::ArgumentsParser(argc, argv)
                        .addArgument("input,i", "input card", &input_card)
                        .addOptionalArgument("point,p", "point to test", &user_point, point_t(12, 0.3))
                        .addOptionalArgument("scan,s", "number of values to scan for a non-zero ME", &num_scans, 0)
                        .parse();

  if (args.debugging())
    CG_LOG_LEVEL(debugInsideLoop);

  cepgen::Generator gen;
  gen.parseRunParameters(input_card);
  gen.runParameters().process().initialise();
  CG_DEBUG("main") << gen.runParameters();

  const auto ndim = gen.runParameters().process().ndim();
  vector<point_t> points;
  if (num_scans > 0) {
    for (const auto& range : cepgen::Limits{0., 1.}.generate(num_scans))
      points.emplace_back(ndim, range);
  } else
    points = vector<point_t>{user_point};
  CG_DEBUG("main") << points;

  vector<pair<point_t, double> > values;
  for (auto& point : points) {
    if (point.size() < 2)
      point = point_t(ndim, point[0]);
    else if (point.size() != ndim)
      point.resize(ndim);
    const double weight = gen.computePoint(point);
    CG_DEBUG("main") << "point " << point << ": weight=" << weight;
    if (weight > 0.)
      values.emplace_back(make_pair(point, weight));
  }
  CG_LOG << "Points with non-zero values: " << values;

  return 0;
}
