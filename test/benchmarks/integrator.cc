#include <nanobench.h>

#include "CepGen/Generator.h"
#include "CepGen/Integration/FunctionalIntegrand.h"
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Modules/IntegratorFactory.h"

int main() {
  cepgen::initialise();
  cepgen::FunctionalIntegrand integrand("x+y^2+z^3", {"x", "y", "z"}, "ROOT");

  for (const auto& integrator_name : {"Vegas", "MISER", "plain"})
    ankerl::nanobench::Bench().run(integrator_name, [&] {
      auto integr = cepgen::IntegratorFactory::get().build(integrator_name);
      integr->setIntegrand(integrand);
      double result, unc;
      integr->integrate(result, unc);
    });

  return 0;
}
