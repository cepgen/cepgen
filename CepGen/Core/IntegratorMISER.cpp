#include "CepGen/Core/Integrator.h"
#include "CepGen/Core/GridParameters.h"
#include "CepGen/Core/Exception.h"

#include "CepGen/Processes/Process.h"

#include "CepGen/Modules/EventModifier.h"
#include "CepGen/Modules/ExportModule.h"

#include "CepGen/Parameters.h"

#include "CepGen/Event/Event.h"
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
      gsl_monte_miser_params miser_params_;
  };

  IntegratorMISER::IntegratorMISER( const ParametersList& params ) :
    Integrator( params )
  {
    miser_params_.estimate_frac = params.get<double>( "estimateFraction", miser_params_.estimate_frac );
    miser_params_.min_calls = params.get<int>( "minCalls", miser_params_.min_calls );
    miser_params_.min_calls_per_bisection = params.get<int>( "minCallsPerBisection", miser_params_.min_calls_per_bisection );
    miser_params_.alpha = params.get<double>( "alpha", miser_params_.alpha );
    miser_params_.dither = params.get<double>( "dither", miser_params_.dither );

    //--- a bit of printout for debugging

    CG_DEBUG( "Integrator:build" ) << "MISER parameters:\n\t"
      << "Number of calls: " << miser_params_.min_calls << ", "
      << "per bisection: " << miser_params_.min_calls_per_bisection << ",\n\t"
      << "Estimate fraction: " << miser_params_.estimate_frac << ",\n\t"
      << "α-value: " << miser_params_.alpha << ",\n\t"
      << "Dither: " << miser_params_.dither << ".";
  }

  void
  IntegratorMISER::integrate( double& result, double& abserr )
  {
    int res = -1;

    //--- integration bounds
    std::vector<double> x_low( function_->dim, 0. ), x_up( function_->dim, 1. );

    //--- launch integration
    std::unique_ptr<gsl_monte_miser_state,void(*)( gsl_monte_miser_state* )>
      mis_state( gsl_monte_miser_alloc( function_->dim ), gsl_monte_miser_free );
    gsl_monte_miser_params_set( mis_state.get(), &miser_params_ );
    res = gsl_monte_miser_integrate( function_.get(),
      &x_low[0], &x_up[0],
      function_->dim, ncvg_,
      rng_.get(), mis_state.get(),
      &result, &abserr );

    result_ = result;
    err_result_ = abserr;

    for ( auto& mod : input_params_->eventModifiersSequence() )
      mod->setCrossSection( result, abserr );
    for ( auto& mod : input_params_->outputModulesSequence() )
      mod->setCrossSection( result, abserr );

    if ( res != GSL_SUCCESS )
      throw CG_FATAL( "Integrator:integrate" )
        << "Error while performing the integration!\n\t"
        << "GSL error: " << gsl_strerror( res ) << ".";
  }
}

