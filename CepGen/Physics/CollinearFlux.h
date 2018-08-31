#ifndef CepGen_Physics_CollinearFlux_h
#define CepGen_Physics_CollinearFlux_h

#include <gsl/gsl_integration.h>
#include <memory>

namespace CepGen
{
  class Limits;
  class StructureFunctions;
  class HeavyIon;
  enum class KTFlux;
  struct FluxArguments
  {
    double x, mx;
    KTFlux flux_type;
    StructureFunctions* str_functions;
    HeavyIon* heavy_ion;
  };
  class CollinearFlux
  {
    public:
      CollinearFlux( const KTFlux& flux_type, const Limits& range, StructureFunctions* str_fun = nullptr );
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

