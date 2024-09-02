/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024  Laurent Forthomme
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
#include "CepGen/Modules/GeneratorWorkerFactory.h"
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

  int num_epochs, min_epochs_iterations, num_events;
  string filename, integrator_name;
  vector<string> processes, generators, outputs;
  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("epochs,e", "number of epochs to try", &num_epochs, 50)
      .addOptionalArgument("epochs-iterations,I", "minimum epochs iterations", &min_epochs_iterations, 5'000)
      .addOptionalArgument("processes,p", "process to benchmark", &processes, vector<string>{"lpair"})
      .addOptionalArgument(
          "generators,g", "event generators to benchmark", &generators, cepgen::GeneratorWorkerFactory::get().modules())
      .addOptionalArgument("num-events,n", "number of events to generate on benchmark", &num_events, 1000)
      .addOptionalArgument("outputs,o", "output formats (html, csv, json, pyperf)", &outputs, vector<string>{"html"})
      .addOptionalArgument("integrator,i", "integrator to use prior to event generation", &integrator_name, "Vegas")
      .addOptionalArgument("filename,f",
                           "output filename",
                           &filename,
                           fs::path(cepgen::utils::env::get("CEPGEN_PATH", ".")) / "benchmark_generator_process")
      .parse();

  ankerl::nanobench::Bench bench;
  bench.title("CepGen v" + cepgen::version::tag + " (" + cepgen::version::extended + ")")
      .epochs(num_epochs)
      .minEpochIterations(min_epochs_iterations);
  for (const auto& process : processes) {
    bench.context("process", process);
    gen.runParameters().setProcess(cepgen::ProcessFactory::get().build(process));
    gen.runParameters().process().kinematics().setParameters(cepgen::ParametersList()
                                                                 .set<vector<int> >("pdgIds", {2212, 2212})
                                                                 .set<double>("sqrtS", 13.6e3)
                                                                 .set<int>("mode", 1)
                                                                 .set<double>("ptmin", 25.));
    gen.setIntegrator(cepgen::IntegratorFactory::get().build(integrator_name));
    for (const auto& generator_name : generators) {
      gen.integrate();  // prepare the grid first
      bench.context("generator", generator_name)
          .run(process + "+" + generator_name, [&gen, &generator_name, &num_events] {
            gen.runParameters().generation().setParameters(cepgen::ParametersList().set(
                "worker", cepgen::GeneratorWorkerFactory::get().describeParameters(generator_name).parameters()));
            gen.generate(num_events);
          });
    }
  }
  render_benchmark(bench, filename, outputs);

  return 0;
}
