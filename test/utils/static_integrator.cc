#include <functional>

#include "CepGen/Generator.h"
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Test.h"

using namespace std;

int main(int argc, char* argv[]) {
  cepgen::initialise();

  string integrator_name;

  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("integrator,i", "integration algorithm", &integrator_name, "Vegas")
      .parse();

  auto integrator_parameters = cepgen::ParametersList().setName(integrator_name);
  {
    constexpr auto alpha = 1.;
    CG_TEST_EQUIV(
        cepgen::Integrator::integrate(
            [&alpha](const vector<double>& vars) -> double { return log(alpha * vars.at(0)) / sqrt(vars.at(0)); },
            integrator_parameters,
            1),
        -4.,
        "standard 1D integration");
  }
  {
    CG_TEST_SET_PRECISION(0.01);
    CG_TEST_EQUIV(cepgen::Integrator::integrate(
                      [](const vector<double>& vars) -> double {
                        double A = 1. / (M_PI * M_PI * M_PI);
                        return A / (1. - cos(vars.at(0)) * cos(vars.at(1)) * cos(vars.at(2)));
                      },
                      integrator_parameters,
                      {cepgen::Limits(0., M_PI), cepgen::Limits(0., M_PI), cepgen::Limits(0., M_PI)}),
                  1.3932039296856768591842462603255,
                  "standard 3D integration");
    CG_TEST_RESET_PRECISION();
  }
  CG_TEST_SUMMARY;
}
