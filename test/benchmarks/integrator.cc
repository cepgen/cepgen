#include <benchmark/benchmark.h>

#include "CepGen/Generator.h"
#include "CepGen/Integration/FunctionalIntegrand.h"
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Modules/IntegratorFactory.h"

cepgen::Integrand* integrand = nullptr;

void run_integrator(benchmark::State& state, const std::string& integr_name) {
  auto integr = cepgen::IntegratorFactory::get().build(integr_name);
  integr->setIntegrand(*integrand);
  for (auto st : state) {
    double result, unc;
    integr->integrate(result, unc);
  }
}

void run_integrator_vegas(benchmark::State& state) { run_integrator(state, "Vegas"); }
void run_integrator_plain(benchmark::State& state) { run_integrator(state, "plain"); }
void run_integrator_miser(benchmark::State& state) { run_integrator(state, "MISER"); }

BENCHMARK(run_integrator_vegas);
BENCHMARK(run_integrator_plain);
BENCHMARK(run_integrator_miser);

//BENCHMARK_MAIN();
int main(int argc, char* argv[]) {
  cepgen::initialise();
  integrand = new cepgen::FunctionalIntegrand("x+y^2+z^3", {"x", "y", "z"}, "ROOT");
  ::benchmark::Initialize(&argc, argv);
  ::benchmark::RunSpecifiedBenchmarks();
}
