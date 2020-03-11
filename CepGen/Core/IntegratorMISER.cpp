#include "CepGen/Core/Integrator.h"
#include "CepGen/Core/GridParameters.h"
#include "CepGen/Core/Exception.h"

#include "CepGen/Modules/IntegratorFactory.h"

#include "CepGen/Parameters.h"
#include "CepGen/Utils/String.h"

#include <gsl/gsl_monte_miser.h>

namespace cepgen
{
  /// MISER integration algorithm developed by W.H. Press and G.R. Farrar, as documented in \cite Press:1989vk.
  class IntegratorMISER : public Integrator
  {
    public:
      IntegratorMISER( const ParametersList& );

      void integrate( double&, double& ) override;

    private:
      std::unique_ptr<gsl_monte_miser_state,void(*)( gsl_monte_miser_state* )> miser_state_;
  };

  IntegratorMISER::IntegratorMISER( const ParametersList& params ) :
    Integrator( params ),
    miser_state_( gsl_monte_miser_alloc( function_->dim ), gsl_monte_miser_free )
  {
    gsl_monte_miser_params miser_params;
    miser_params.estimate_frac = params.get<double>( "estimateFraction", miser_params.estimate_frac );
    miser_params.min_calls = params.get<int>( "minCalls", miser_params.min_calls );
    miser_params.min_calls_per_bisection = params.get<int>( "minCallsPerBisection", miser_params.min_calls_per_bisection );
    miser_params.alpha = params.get<double>( "alpha", miser_params.alpha );
    miser_params.dither = params.get<double>( "dither", miser_params.dither );
    gsl_monte_miser_params_set( miser_state_.get(), &miser_params );

    //--- a bit of printout for debugging

    CG_DEBUG( "Integrator:build" ) << "MISER parameters:\n\t"
      << "Number of calls: " << miser_params.min_calls << ", "
      << "per bisection: " << miser_params.min_calls_per_bisection << ",\n\t"
      << "Estimate fraction: " << miser_params.estimate_frac << ",\n\t"
      << "α-value: " << miser_params.alpha << ",\n\t"
      << "Dither: " << miser_params.dither << ".";
  }

  void
  IntegratorMISER::integrate( double& result, double& abserr )
  {
    //--- integration bounds
    std::vector<double> x_low( function_->dim, 0. ), x_up( function_->dim, 1. );

    //--- launch integration
    int res = gsl_monte_miser_integrate( function_.get(),
      &x_low[0], &x_up[0],
      function_->dim, ncvg_,
      rng_.get(), miser_state_.get(),
      &result, &abserr );

    if ( res != GSL_SUCCESS )
      throw CG_FATAL( "Integrator:integrate" )
        << "Error while performing the integration!\n\t"
        << "GSL error: " << gsl_strerror( res ) << ".";

    result_ = result;
    err_result_ = abserr;
  }
}

REGISTER_INTEGRATOR( "MISER", IntegratorMISER )
