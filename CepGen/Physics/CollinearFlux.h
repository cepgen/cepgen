#ifndef CepGen_Physics_CollinearFlux_h
#define CepGen_Physics_CollinearFlux_h

#include <gsl/gsl_integration.h>

#include <memory>

namespace cepgen {
  class Limits;
  namespace formfac {
    class Parameterisation;
  }
  class HeavyIon;
  enum class KTFlux;
  struct FluxArguments {
    double x, mi2, mf2;
    KTFlux flux_type;
    formfac::Parameterisation* form_factors;
    HeavyIon* heavy_ion;
  };
  class CollinearFlux {
  public:
    CollinearFlux(formfac::Parameterisation* form_fac, const Limits& kt2_range);
    CollinearFlux(HeavyIon* hi, const Limits& kt2_range);
    double operator()(double x, double mx, const KTFlux& flux_type) const;

  private:
    struct gsl_integration_fixed_workspace_del {
      void operator()(gsl_integration_fixed_workspace* int_wsp) { gsl_integration_fixed_free(int_wsp); }
    };
    std::unique_ptr<gsl_integration_fixed_workspace, gsl_integration_fixed_workspace_del> workspace_;
    std::unique_ptr<FluxArguments> params_;
    mutable gsl_function function_;
  };
}  // namespace cepgen

#endif
