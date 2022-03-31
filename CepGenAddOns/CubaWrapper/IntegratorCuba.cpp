#include "CepGen/Core/Exception.h"
#include "CepGen/Integration/Integrand.h"
#include "CepGenAddOns/CubaWrapper/IntegratorCuba.h"

namespace cepgen {
  IntegratorCuba::IntegratorCuba(const ParametersList& params)
      : Integrator(params),
        ncomp_(steer<int>("ncomp")),
        nvec_(steer<int>("nvec")),
        epsrel_(steer<double>("epsrel")),
        epsabs_(steer<double>("epsabs")),
        mineval_(steer<int>("mineval")),
        maxeval_(steer<int>("maxeval")),
        verbose_(steer<int>("verbose")) {}

  void IntegratorCuba::setIntegrand(Integrand& integr) {
    Integrator::setIntegrand(integr);
    gIntegrand = integrand_;
  }

  ParametersDescription IntegratorCuba::description() {
    auto desc = Integrator::description();
    desc.setDescription("Cuba generic integration algorithm");
    desc.add<int>("ncomp", 1).setDescription("number of components of the integrand");
    desc.add<int>("nvec", 1).setDescription("number of samples received by the integrand");
    desc.add<double>("epsrel", 1.e-3).setDescription("requested relative accuracy");
    desc.add<double>("epsabs", 1.e-12).setDescription("requested absolute accuracy");
    desc.add<int>("mineval", 0).setDescription("minimum number of integrand evaluations required");
    desc.add<int>("maxeval", 50000).setDescription("(approximate) maximum number of integrand evaluations allowed");
    desc.add<int>("verbose", 0);
    return desc;
  }

  int cuba_integrand(const int* ndim, const double xx[], const int* /*ncomp*/, double ff[], void* /*userdata*/) {
    if (!gIntegrand)
      throw CG_FATAL("cuba_integrand") << "Integrand not set for the Cuba algorithm!";
    ff[0] = gIntegrand->eval(std::vector<double>(xx, xx + *ndim));
    return 0;
  }
}  // namespace cepgen
