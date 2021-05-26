#ifndef CepGenAddOns_CubaWrapper_IntegratorCuba_h
#define CepGenAddOns_CubaWrapper_IntegratorCuba_h

#include "CepGen/Integration/Integrator.h"
#include "CepGen/Integration/Integrand.h"

namespace cepgen {
  /// Cuba integration algorithm
  class IntegratorCuba : public Integrator {
  public:
    explicit IntegratorCuba(const ParametersList& params)
        : Integrator(params),
          nvec_(params.get<int>("NVEC", 1)),
          epsrel_(params.get<double>("EPSREL", 1.e-3)),
          epsabs_(params.get<double>("EPSABS", 1.e-12)),
          mineval_(params.get<int>("MINEVAL", 0)),
          maxeval_(params.get<int>("MAXEVAL", 50000)),
          verbose_(params.get<int>("verbose", 1)) {}
    static std::string description() { return "Cuba generic integration algorithm"; }

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
