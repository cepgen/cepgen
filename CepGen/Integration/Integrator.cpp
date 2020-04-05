#include "CepGen/Integration/Integrator.h"
#include "CepGen/Integration/Integrand.h"

#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/EventModifier.h"
#include "CepGen/Core/ExportModule.h"
#include "CepGen/Core/Exception.h"

#include "CepGen/Parameters.h"
#include "CepGen/Processes/Process.h"
#include "CepGen/Event/Event.h"

#include "CepGen/Utils/String.h"
#include "CepGen/Utils/ProgressBar.h"

namespace cepgen
{
  Integrator::Integrator( const ParametersList& params ) :
    params_( params ),
    name_( params.name<std::string>() ),
    seed_( params.get<int>( "seed", time( nullptr ) ) ),
    verbosity_( params.get<int>( "verbose", 1 ) ),
    initialised_( false )
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
      << "Random numbers generator: " << gsl_rng_name( rng_.get() ) << ".\n\t"
      << "Seed: " << seed_ << ".";
  }

  void
  Integrator::setIntegrand( Integrand& integr )
  {
    integrand_ = &integr;
    //--- specify the integrand through the GSL wrapper
    auto func = [=]( double* x, size_t ndim, void* ) -> double {
      return integrand_->eval( std::vector<double>( x, x+ndim ) );
    };
    function_.reset( new gsl_monte_function_wrapper<decltype( func )>( func, integrand_->size() ) );

    CG_DEBUG( "Integrator:integrand" )
      << "Number of integration dimensions: " << function_->dim << ".";

    //--- force the reinitialisation
    initialised_ = false;
  }

  //------------------------------------------------------------------------------------------------
  // helper / alias methods
  //------------------------------------------------------------------------------------------------

  size_t
  Integrator::size() const
  {
    if ( !function_ )
      throw CG_FATAL( "Integrator:size" )
        << "Trying to retrieve phase space size on an unitialised integrand!";
    return function_->dim;
  }

  double
  Integrator::eval( const std::vector<double>& x ) const
  {
    if ( !function_ )
      throw CG_FATAL( "Integrator:eval" )
        << "Trying to evaluate the weight on a phase space point "
        << "on an unitialised integrand!";
    return integrand_->eval( x );
  }

  double
  Integrator::uniform() const
  {
    if ( !rng_ )
      throw CG_FATAL( "Integrator:uniform" )
        << "Random number generator has not been initialised!";
    return gsl_rng_uniform( rng_.get() );
  }
}

