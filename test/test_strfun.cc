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

#include "CepGen/Core/Exception.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Utils/ArgumentsParser.h"

using namespace std;

int main(int argc, char* argv[]) {
  vector<int> strfun_type;
  int num_points;
  vector<double> q2in, xbjin;
  string output_file;

  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("sf,s", "structure functions modelling", &strfun_type, vector<int>{301})
      .addOptionalArgument("q2,q", "parton virtuality (GeV^2)", &q2in, vector<double>{})
      .addOptionalArgument("xbj,x", "Bjorken-x", &xbjin, vector<double>{})
      .addOptionalArgument("npoints,n", "number of x-points to scan", &num_points, 100)
      .addOptionalArgument("output,o", "output file name", &output_file, "flux.scan.output.txt")
      .parse();

  cepgen::initialise();
  vector<double> q2vals, xbjvals;

  ofstream out(output_file);
  out << "# structure functions:\n";
  vector<std::unique_ptr<cepgen::strfun::Parameterisation> > params;
  for (const auto& type : strfun_type) {
    params.emplace_back(cepgen::strfun::StructureFunctionsFactory::get().build(type));
    auto& sf = *params.rbegin();
    out << "# * " << sf.get() << "\n";
  }
  if (q2in.empty())
    throw CG_FATAL("main") << "At least one value of Q^2 is required!";
  else if (q2in.size() == 2) {  // min-max
    q2vals.clear();
    for (int i = 0; i <= num_points; ++i)
      q2vals.emplace_back(q2in[0] + i * (q2in[1] - q2in[0]) / num_points);
  } else
    q2vals = q2in;

  if (xbjin.empty())
    throw CG_FATAL("main") << "At least one value of x_Bj is required!";
  else if (xbjin.size() == 2) {  // min-max
    xbjvals.clear();
    for (int i = 0; i <= num_points; ++i)
      xbjvals.emplace_back(xbjin[0] + i * (xbjin[1] - xbjin[0]) / num_points);
  } else
    xbjvals = xbjin;

  out << "# q2\txbj\tF_2\tF_L\n";

  for (const auto& xbj : xbjvals)
    for (const auto& q2 : q2vals) {
      out << q2 << "\t" << xbj;
      for (auto& sf : params) {
        auto& sfval = (*sf)(xbj, q2);
        sfval.computeFL(xbj, q2);
        out << "\t" << sfval.F2 << "\t" << sfval.FL;
      }
      out << "\n";
    }

  CG_LOG << "Scan written in \"" << output_file << "\".";
  out.close();

  return 0;
}
