#include "CepGen/Physics/CollinearFlux.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/StructureFunctions/StructureFunctions.h"

#include <gsl/gsl_errno.h>

namespace CepGen
{
  double
  unintegrated_flux( double kt2, void* params )
  {
    const auto args = (FluxArguments*)params;
    return Process::GenericKTProcess::flux( args->flux_type, args->x, kt2, *args->str_functions, args->mx );
  }

  CollinearFlux::CollinearFlux( const Process::GenericKTProcess::Flux& flux_type, const Limits& range, StructureFunctions* str_fun ) :
    workspace_( gsl_integration_fixed_alloc( gsl_integration_fixed_legendre, 100, range.min(), range.max(), 0., 0. ) ),
    params_( new FluxArguments{ flux_type, str_fun, 0., 0. } ),
    function_( { &unintegrated_flux, (void*)params_.get() } )
  {}

  double
  CollinearFlux::operator()( double x, double mx ) const
  {
    double result = 0.;
    params_->x = x;
    params_->mx = mx;
    const int res = gsl_integration_fixed( &function_, &result, workspace_.get() );
    if ( res != GSL_SUCCESS )
      CG_ERROR( "CollinearFlux" ) << gsl_strerror( res );
    return result;
  }
}

