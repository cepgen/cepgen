#ifndef CepGen_Physics_CollinearFlux_h
#define CepGen_Physics_CollinearFlux_h

#include "CepGen/Processes/GenericKTProcess.h"

#include <gsl/gsl_integration.h>
#include <memory>

namespace CepGen
{
  class Limits;
  struct FluxArguments
  {
    Process::GenericKTProcess::Flux flux_type;
    StructureFunctions* str_functions;
    double x, mx;
  };
  class CollinearFlux
  {
    public:
      CollinearFlux( const Process::GenericKTProcess::Flux& flux_type, const Limits& range, StructureFunctions* str_fun = nullptr );
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

