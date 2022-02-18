#ifndef CepGenAddOns_CubaWrapper_IntegratorCuba_h
#define CepGenAddOns_CubaWrapper_IntegratorCuba_h

#include "CepGen/Integration/Integrand.h"
#include "CepGen/Integration/Integrator.h"

namespace cepgen {
  /// Cuba integration algorithm
  class IntegratorCuba : public Integrator {
  public:
    explicit IntegratorCuba(const ParametersList& params)
        : Integrator(params),
          nvec_(steer<int>("nvec")),
          epsrel_(steer<double>("epsrel")),
          epsabs_(steer<double>("epsabs")),
          mineval_(steer<int>("mineval")),
          maxeval_(steer<int>("maxeval")),
          verbose_(steer<int>("verbose")) {}

    static ParametersDescription description() {
      auto desc = Integrator::description();
      desc.setDescription("Cuba generic integration algorithm");
      desc.add<int>("nvec", 1).setDescription("number of samples received by the integrand");
      desc.add<double>("epsrel", 1.e-3).setDescription("requested relative accuracy");
      desc.add<double>("epsabs", 1.e-12).setDescription("requested absolute accuracy");
      desc.add<int>("mineval", 0).setDescription("minimum number of integrand evaluations required");
      desc.add<int>("maxeval", 50000).setDescription("(approximate) maximum number of integrand evaluations allowed");
      desc.add<int>("verbose", 1);
      return desc;
    }

  protected:
    int nvec_;
    double epsrel_, epsabs_;
    int mineval_, maxeval_;
    int verbose_;
  };

  static Integrand* gIntegrand = nullptr;
  inline int cuba_integrand(const int* ndim, const double xx[], const int* /*ncomp*/, double ff[], void* /*userdata*/) {
    ff[0] = gIntegrand->eval(std::vector<double>(xx, xx + *ndim));
    return 0;
  }
}  // namespace cepgen

#endif
