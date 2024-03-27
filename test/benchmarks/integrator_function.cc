/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023  Laurent Forthomme
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

#include "CepGen/Generator.h"
#include "CepGen/Integration/FunctionalIntegrand.h"
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Modules/FunctionalFactory.h"
#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Environment.h"
#include "CepGen/Utils/Filesystem.h"
#include "CepGen/Utils/Logger.h"
#include "CepGen/Version.h"
#include "nanobench_interface.h"

using namespace std;

int main(int argc, char* argv[]) {
  cepgen::initialise();

  int num_epochs;
  string filename;
  vector<string> functional_parsers, integrators, outputs;
  bool python_integ;
  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("epochs,e", "number of epochs to try", &num_epochs, 10)
      .addOptionalArgument("functionals,F",
                           "functional parsers to benchmark",
                           &functional_parsers,
                           cepgen::FunctionalFactory::get().modules())
      .addOptionalArgument(
          "integrators,i", "integrators to benchmark", &integrators, cepgen::IntegratorFactory::get().modules())
      .addOptionalArgument("outputs,o", "output formats (html, csv, json, pyperf)", &outputs, vector<string>{"html"})
      .addOptionalArgument("filename,f",
                           "output filename",
                           &filename,
                           fs::path(cepgen::utils::env::get("CEPGEN_PATH", ".")) / "benchmark_integrator_function")
      .addOptionalArgument("python,p", "also add python integrator?", &python_integ, false)
      .parse();

  ankerl::nanobench::Bench bench;
  bench.title("CepGen v" + cepgen::version::tag + " (" + cepgen::version::extended + ")").epochs(num_epochs);
  for (const auto& functional_parser : functional_parsers) {
    bench.context("functional", functional_parser);
    cepgen::FunctionalIntegrand integrand("x+y^2+z^3", {"x", "y", "z"}, functional_parser);
    for (const auto& integrator_name : integrators) {
      if (integrator_name == "python" && !python_integ)  // skip the python integrators test unless required
        continue;
      bench.context("integrator", integrator_name).run(functional_parser + "+" + integrator_name, [&] {
        auto integr = cepgen::IntegratorFactory::get().build(integrator_name);
        integr->integrate(integrand);
      });
    }
  }
  render_benchmark(bench, filename, outputs);

  return 0;
}
