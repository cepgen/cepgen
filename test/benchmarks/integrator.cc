#include <nanobench.h>

#include <fstream>

#include "CepGen/Generator.h"
#include "CepGen/Integration/FunctionalIntegrand.h"
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Modules/FunctionalFactory.h"
#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Version.h"

using namespace std;

int main(int argc, char* argv[]) {
  cepgen::initialise();

  int num_epochs;
  vector<string> functional_parsers, integrators;
  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("epochs,e", "number of epochs to try", &num_epochs, 20)
      .addOptionalArgument("functionals,f",
                           "functional parsers to benchmark",
                           &functional_parsers,
                           cepgen::utils::FunctionalFactory::get().modules())
      .addOptionalArgument(
          "integrators,i", "integrators to benchmark", &integrators, cepgen::IntegratorFactory::get().modules())
      .parse();

  ofstream out_file("benchmark.html");

  ankerl::nanobench::Bench bench;
  bench.title(("CepGen v" + cepgen::version::tag + " (" + cepgen::version::extended + ")").data()).epochs(num_epochs);

  for (const auto& functional_parser : functional_parsers) {
    bench.context("functional", functional_parser.data());
    cepgen::FunctionalIntegrand integrand("x+y^2+z^3", {"x", "y", "z"}, functional_parser);
    for (const auto& integrator_name : integrators)
      bench.context("integrator", integrator_name.data()).run((functional_parser + "+" + integrator_name).data(), [&] {
        auto integr = cepgen::IntegratorFactory::get().build(integrator_name);
        integr->setIntegrand(integrand);
        double result, unc;
        integr->integrate(result, unc);
      });
    bench.render(ankerl::nanobench::templates::htmlBoxplot(), out_file);
  }

  return 0;
}
