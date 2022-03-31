#ifndef CepGenAddOns_CubaWrapper_IntegratorCuba_h
#define CepGenAddOns_CubaWrapper_IntegratorCuba_h

#include "CepGen/Integration/Integrator.h"

namespace cepgen {
  class Integrand;
  /// Cuba integration algorithm
  class IntegratorCuba : public Integrator {
  public:
    explicit IntegratorCuba(const ParametersList&);

    static ParametersDescription description();
    void setIntegrand(Integrand&) override;

  protected:
    int ncomp_, nvec_;
    double epsrel_, epsabs_;
    int mineval_, maxeval_;
    int verbose_;
  };

  int cuba_integrand(const int* ndim, const double xx[], const int* /*ncomp*/, double ff[], void* /*userdata*/);
}  // namespace cepgen

#endif
