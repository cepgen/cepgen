#include "CepGen/Integration/Integrator.h"

#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/Exception.h"

#include "CepGen/Parameters.h"
#include "CepGen/Processes/Process.h"
#include "CepGen/Event/Event.h"

#include "CepGen/Modules/EventModifier.h"
#include "CepGen/Modules/ExportModule.h"

#include "CepGen/Utils/String.h"
#include "CepGen/Utils/ProgressBar.h"

namespace cepgen
{
  Integrator::Integrator( const ParametersList& params ) :
    params_( params ),
    name_( params.name<std::string>() ),
    seed_( params.get<int>( "seed", time( nullptr ) ) ),
    initialised_( false ), ps_bin_( INVALID_BIN )
  {
    //--- initialise the random number generator
    const auto& rng_type = params.get<int>( "rngEngine" );
    gsl_rng_type* rng_engine = nullptr;
    switch ( rng_type ) {
      case 0: default: rng_engine = (gsl_rng_type*)gsl_rng_mt19937; break;
      case 1: rng_engine = (gsl_rng_type*)gsl_rng_taus2; break;
      case 2: rng_engine = (gsl_rng_type*)gsl_rng_gfsr4; break;
      case 3: rng_engine = (gsl_rng_type*)gsl_rng_ranlxs0; break;
    }
    if ( !rng_engine  )
      throw CG_FATAL( "Integrator:build" )
        << "Random number generator engine not set!";

    rng_.reset( gsl_rng_alloc( rng_engine ) );
    gsl_rng_set( rng_.get(), seed_ );

    //--- a bit of printout for debugging

    CG_DEBUG( "Integrator:build" )
      << "Number of function calls: " << ncvg_ << ",\n\t"
      << "Random numbers generator: " << gsl_rng_name( rng_.get() ) << ".";
  }

  void
  Integrator::setFunction( unsigned int ndim, double integrand( double*, size_t, void* ), Parameters& params )
  {
    input_params_ = &params;
    function_.reset( new gsl_monte_function{ integrand, ndim, (void*)input_params_ } );

    CG_DEBUG( "Integrator:function" )
      << "Number of integration dimensions: " << function_->dim << ".";

    initialised_ = false;
  }

  //------------------------------------------------------------------------------------------------
  // helper / alias methods
  //------------------------------------------------------------------------------------------------

  size_t
  Integrator::size() const
  {
    if ( !function_ )
      return 0;
    return function_->dim;
  }

  double
  Integrator::eval( const std::vector<double>& x )
  {
    //--- by default, no grid treatment
    return function_->f( (double*)&x[0], function_->dim, (void*)input_params_ );
  }

  double
  Integrator::uniform() const
  {
    return gsl_rng_uniform( rng_.get() );
  }
}

