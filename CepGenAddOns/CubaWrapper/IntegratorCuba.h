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
          nvec_(steer<int>("NVEC")),
          epsrel_(steer<double>("EPSREL")),
          epsabs_(steer<double>("EPSABS")),
          mineval_(steer<int>("MINEVAL")),
          maxeval_(steer<int>("MAXEVAL")),
          verbose_(steer<int>("verbose")) {}

    static ParametersDescription description() {
      auto desc = Integrator::description();
      desc.setDescription("Cuba generic integration algorithm");
      desc.add<int>("NVEC", 1);
      desc.add<double>("EPSREL", 1.e-3);
      desc.add<double>("EPSABS", 1.e-12);
      desc.add<int>("MINEVAL", 0);
      desc.add<int>("MAXEVAL", 50000);
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
