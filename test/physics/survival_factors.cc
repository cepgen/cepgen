/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2026  Laurent Forthomme
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
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Test.h"

using namespace std;

int main(int argc, char* argv[]) {
  double num_sigma, s2_el;
  string proc_name;

  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("process,p", "process to compute", &proc_name, "lpair")
      .addOptionalArgument("num-sigma,n", "max. number of std.dev.", &num_sigma, 3.)
      .addOptionalArgument("s2-el,e", "elastic-elastic survival factor value", &s2_el, 0.9)
      .parse();

  cepgen::Generator gen;
  auto compute_xsec = [&](int mode, const cepgen::ParametersList& survival) {
    auto process_parameters = cepgen::ParametersList().set<int>("pair", 13);
    if (!survival.empty())
      process_parameters.set("survivalFactor", survival);
    gen.runParameters().setProcess(cepgen::ProcessFactory::get().build(proc_name, process_parameters));
    gen.runParameters().process().kinematics().setParameters(
        cepgen::ParametersList()
            .set<double>("sqrtS", 13.e3)
            .set<double>("ptmin", 5.)
            .set<int>("mode", mode)
            .set<cepgen::ParametersList>(
                "structureFunctions",
                cepgen::StructureFunctionsFactory::get().describeParameters("SuriYennie").parameters())
            .set<cepgen::Limits>("eta", {-2.5, 2.5})
            .set<cepgen::Limits>("mx", {1.07, 1000.}));
    return gen.computeXsection();
  };
  {  // elastic-elastic
    const auto cs_none = compute_xsec(1, cepgen::ParametersList());
    const auto cs_constant = compute_xsec(1, cepgen::ParametersList().setName("constant").set<double>("value", s2_el));
    CG_TEST_VALUES(cs_none * s2_el, cs_constant, num_sigma, "el-el: constant S^2 rescales the cross section");
  }

  CG_TEST_SUMMARY;
}
