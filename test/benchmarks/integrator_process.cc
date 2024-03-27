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
#include "CepGen/Utils/Environment.h"
#include "CepGen/Utils/Filesystem.h"
#include "CepGen/Version.h"
#include "nanobench_interface.h"

using namespace std;

int main(int argc, char* argv[]) {
  cepgen::Generator gen;

  int num_epochs;
  string filename;
  vector<string> processes, integrators, outputs;
  bool python_integ;
  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("epochs,e", "number of epochs to try", &num_epochs, 5)
      .addOptionalArgument("processes,p", "process to benchmark", &processes, vector<string>{"lpair"})
      .addOptionalArgument(
          "integrators,i", "integrators to benchmark", &integrators, cepgen::IntegratorFactory::get().modules())
      .addOptionalArgument("outputs,o", "output formats (html, csv, json, pyperf)", &outputs, vector<string>{"html"})
      .addOptionalArgument("filename,f",
                           "output filename",
                           &filename,
                           fs::path(cepgen::utils::env::get("CEPGEN_PATH", ".")) / "benchmark_integrator_process")
      .addOptionalArgument("python,p", "also add python integrator?", &python_integ, false)
      .parse();

  ankerl::nanobench::Bench bench;
  bench.title("CepGen v" + cepgen::version::tag + " (" + cepgen::version::extended + ")").epochs(num_epochs);
  for (const auto& process : processes) {
    bench.context("process", process);
    gen.runParameters().setProcess(cepgen::ProcessFactory::get().build(process));
    gen.runParameters().process().kinematics().setParameters(cepgen::ParametersList()
                                                                 .set<vector<int> >("pdgIds", {2212, 2212})
                                                                 .set<double>("sqrtS", 13.6e3)
                                                                 .set<int>("mode", 1)
                                                                 .set<double>("ptmin", 25.));
    for (const auto& integrator_name : integrators) {
      if (integrator_name == "python" && !python_integ)  // skip the python integrators test unless required
        continue;
      bench.context("integrator", integrator_name).run(process + "+" + integrator_name, [&] {
        gen.setIntegrator(cepgen::IntegratorFactory::get().build(integrator_name));
        gen.computeXsection();
      });
    }
  }
  render_benchmark(bench, filename, outputs);

  return 0;
}
