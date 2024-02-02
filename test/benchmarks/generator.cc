/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2024  Laurent Forthomme
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
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Version.h"
#include "nanobench_interface.h"

using namespace std;

int main(int argc, char* argv[]) {
  cepgen::Generator gen;

  int num_epochs;
  string process;
  vector<string> integrators, outputs;
  bool python_integ;
  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("epochs,e", "number of epochs to try", &num_epochs, 5)
      .addOptionalArgument("process,p", "process to benchmark", &process, "lpair")
      .addOptionalArgument(
          "integrators,i", "integrators to benchmark", &integrators, cepgen::IntegratorFactory::get().modules())
      .addOptionalArgument("outputs,o", "output formats (html, csv, json, pyperf)", &outputs, vector<string>{"html"})
      .addOptionalArgument("python,p", "also add python integrator?", &python_integ, false)
      .parse();

  ankerl::nanobench::Bench bench;
  bench.title("CepGen v" + cepgen::version::tag + " (" + cepgen::version::extended + ")")
      .epochs(num_epochs)
      .context("process", process);

  gen.runParameters().setProcess(cepgen::ProcessFactory::get().build(process));
  auto& kin = gen.runParameters().process().kinematics();
  kin.incomingBeams().positive().setPdgId(2212);
  kin.incomingBeams().negative().setPdgId(2212);
  kin.incomingBeams().setSqrtS(13.e3);
  kin.cuts().central.pt_single.min() = 15.;
  kin.cuts().central.eta_single = {-2.5, 2.5};
  for (const auto& integrator_name : integrators) {
    if (integrator_name == "python" && !python_integ)  // skip the python integrators test unless required
      continue;
    bench.context("integrator", integrator_name).run(process + "+" + integrator_name, [&] {
      gen.setIntegrator(cepgen::IntegratorFactory::get().build(integrator_name));
      gen.computeXsection();
    });
  }
  render_benchmark(bench, outputs);

  return 0;
}
