#ifndef CepGen_Physics_CollinearFlux_h
#define CepGen_Physics_CollinearFlux_h

#include <gsl/gsl_integration.h>
#include <memory>

namespace CepGen
{
  class Limits;
  namespace formfac{ class Parameterisation; }
  class HeavyIon;
  enum class KTFlux;
  struct FluxArguments
  {
    double x, mi2, mf2;
    KTFlux flux_type;
    formfac::Parameterisation* form_factors;
    HeavyIon* heavy_ion;
  };
  class CollinearFlux
  {
    public:
      CollinearFlux( const KTFlux& flux_type, const Limits& range, formfac::Parameterisation* form_fac = nullptr );
      CollinearFlux( const KTFlux& flux_type, const Limits& range, HeavyIon* hi );
      double operator()( double x, double mx ) const;

    private:
      struct gsl_integration_fixed_workspace_del
      {
        void operator()( gsl_integration_fixed_workspace* int_wsp ) {
          gsl_integration_fixed_free( int_wsp );
        }
      };
      std::unique_ptr<gsl_integration_fixed_workspace,gsl_integration_fixed_workspace_del> workspace_;
      std::unique_ptr<FluxArguments> params_;
      mutable gsl_function function_;
  };
}

#endif

