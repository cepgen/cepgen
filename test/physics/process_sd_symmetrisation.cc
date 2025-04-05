/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2025  Laurent Forthomme
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
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/AbortHandler.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Test.h"
#include "CepGen/Utils/Timer.h"

using namespace std;

int main(int argc, char* argv[]) {
  double num_sigma;
  string str_fun, proc_name, integrator;

  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("process,p", "process to compute", &proc_name, "lpair")
      .addOptionalArgument("num-sigma,n", "max. number of std.dev.", &num_sigma, 3.)
      .addOptionalArgument("str-fun,s", "struct.functions modelling", &str_fun, "SuriYennie")
      .addOptionalArgument("integrator,i", "type of integrator used", &integrator, "Vegas")
      .parse();

  [[maybe_unused]] cepgen::utils::Timer timer;
  cepgen::Generator gen;
  gen.runParameters().integrator().setName(integrator);

  [[maybe_unused]] cepgen::utils::AbortHandler abort_handler;

  auto kinematics_parameters =
      cepgen::ParametersList()
          .set<double>("sqrtS", 13.e3)
          .set<cepgen::ParametersList>(
              "structureFunctions", cepgen::StructureFunctionsFactory::get().describeParameters(str_fun).parameters())
          .set<double>("ptmin", 5.)
          .set<int>("mode", 2)  // elastic-inelastic
          .set<cepgen::Limits>("eta", {-2.5, 2.5})
          .set<cepgen::Limits>("mx", {1.07, 1000.});

  cepgen::Value cs_ei_no_symm, cs_ei_symm;
  {
    gen.runParameters().setProcess(
        cepgen::ProcessFactory::get().build(proc_name, cepgen::ParametersList().set<int>("pair", 13)));
    gen.runParameters().process().kinematics().setParameters(kinematics_parameters);
    cs_ei_no_symm = gen.computeXsection();
  }
  {  // inelastic-elastic
    gen.runParameters().setProcess(cepgen::ProcessFactory::get().build(
        proc_name, cepgen::ParametersList().set<int>("pair", 13).set<bool>("symmetrise", true)));
    gen.runParameters().process().kinematics().setParameters(kinematics_parameters);
    cs_ei_symm = gen.computeXsection();
  }
  CG_TEST_VALUES(cs_ei_no_symm * 2., cs_ei_symm, num_sigma, "symmetrised SD == 2*non-symmetrised SD");

  CG_TEST_SUMMARY;
}
