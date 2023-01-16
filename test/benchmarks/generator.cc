#include "CepGen/Generator.h"
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Parameters.h"
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
  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("epochs,e", "number of epochs to try", &num_epochs, 5)
      .addOptionalArgument("process,p", "process to benchmark", &process, "lpair")
      .addOptionalArgument(
          "integrators,i", "integrators to benchmark", &integrators, cepgen::IntegratorFactory::get().modules())
      .addOptionalArgument("outputs,o", "output formats (html, csv, json, pyperf)", &outputs, vector<string>{"html"})
      .parse();

  ankerl::nanobench::Bench bench;
  bench.title("CepGen v" + cepgen::version::tag + " (" + cepgen::version::extended + ")")
      .epochs(num_epochs)
      .context("process", process);

  gen.parametersRef().setProcess(cepgen::proc::ProcessFactory::get().build(process));
  auto& kin = gen.parametersRef().process().kinematics();
  kin.incomingBeams().positive().setPdgId(2212);
  kin.incomingBeams().negative().setPdgId(2212);
  kin.incomingBeams().setSqrtS(13.e3);
  for (const auto& integrator_name : integrators)
    bench.context("integrator", integrator_name).run(process + "+" + integrator_name, [&] {
      gen.setIntegrator(cepgen::IntegratorFactory::get().build(integrator_name));
      double xsec, xsec_unc;
      gen.computeXsection(xsec, xsec_unc);
    });
  render_benchmark(bench, outputs);

  return 0;
}
