#include "CepGen/Core/Integrator.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/GridParameters.h"

#include "CepGen/Modules/IntegratorFactory.h"

#include "CepGen/Parameters.h"
#include "CepGen/Utils/String.h"

#include <gsl/gsl_monte_vegas.h>

#define COORD(s,i,j) ((s)->xi[(i)*(s)->dim + (j)])

namespace cepgen
{
  /// Vegas integration algorithm developed by P. Lepage, as documented in \cite Lepage:1977sw
  class IntegratorVegas : public Integrator
  {
    public:
      IntegratorVegas( const ParametersList& );

      void integrate( double&, double& ) override;

      enum class Mode { importance = 1, importanceOnly = 0, stratified = -1 };

    private:
      void warmup( std::vector<double>&, std::vector<double>&, unsigned int );
      double eval( const std::vector<double>& ) override;

      /// A trivial deleter for the Vegas integrator
      struct gsl_monte_vegas_deleter
      {
        inline void operator()( gsl_monte_vegas_state* state ) { gsl_monte_vegas_free( state ); }
      };
      /// A Vegas integrator state for integration (optional) and/or
      /// "treated" event generation
      std::unique_ptr<gsl_monte_vegas_state,gsl_monte_vegas_deleter> vegas_state_;
      gsl_monte_vegas_params vegas_params_;
      double chisq_cut_;

      /// Collection of integrator parameters
      struct Integration
      {
        gsl_rng_type* rng_engine; ///< Random number generator engine
      };
  };
  std::ostream& operator<<( std::ostream&, const IntegratorVegas::Mode& );

  IntegratorVegas::IntegratorVegas( const ParametersList& params ) :
    Integrator( params ),
    chisq_cut_( params.get<double>( "chiSqCut", 1.5 ) )
  {
    vegas_params_.iterations = params.get<int>( "iterations", 10 );
    vegas_params_.alpha = params.get<double>( "alpha", 1.5 );
    vegas_params_.verbose = params.get<int>( "verbose", 0 );
    vegas_params_.mode = params.get<int>( "mode", (int)Mode::importance );

    //--- output logging
    const auto& log = params.get<std::string>( "loggingOutput", "cerr" );
    if ( log == "cerr" )
      // redirect all debugging information to the error stream
      vegas_params_.ostream = stderr;
    else if ( log == "cout" )
      // redirect all debugging information to the standard stream
      vegas_params_.ostream = stdout;
    else
      vegas_params_.ostream = fopen( log.c_str(), "w" );

    //--- a bit of printout for debugging
    CG_DEBUG( "Integrator:build" ) << "Vegas parameters:\n\t"
      << "Number of iterations in Vegas: " << vegas_params_.iterations << ",\n\t"
      << "α-value: " << vegas_params_.alpha << ",\n\t"
      << "Verbosity: " << vegas_params_.verbose << ",\n\t"
      << "Grid interpolation mode: " << (IntegratorVegas::Mode)vegas_params_.mode << ".";
  }

  void
  IntegratorVegas::integrate( double& result, double& abserr )
  {
    //--- integration bounds
    std::vector<double> x_low( function_->dim, 0. ), x_up( function_->dim, 1. );

    //--- launch integration

    //----- warmup (prepare the grid)
    warmup( x_low, x_up, 25000 );
    //----- integration
    unsigned short it_chisq = 0;
    do {
      int res = gsl_monte_vegas_integrate( function_.get(),
        &x_low[0], &x_up[0],
        function_->dim, 0.2 * ncvg_,
        rng_.get(), vegas_state_.get(),
        &result, &abserr );
      CG_LOG( "Integrator:integrate" )
        << "\t>> at call " << ( ++it_chisq ) << ": "
        << utils::format( "average = %10.6f   "
                          "sigma = %10.6f   chi2 = %4.3f.",
                          result, abserr,
                          gsl_monte_vegas_chisq( vegas_state_.get() ) );
      if ( res != GSL_SUCCESS )
        throw CG_FATAL( "Integrator:integrate" )
          << "Error at iteration #" << it_chisq
          << " while performing the integration!\n\t"
          << "GSL error: " << gsl_strerror( res ) << ".";
    } while ( fabs( gsl_monte_vegas_chisq( vegas_state_.get() )-1. ) > chisq_cut_-1. );
    CG_DEBUG( "Integrator:integrate" )
      << "Vegas grid information:\n\t"
      << "ran for " << vegas_state_->dim << " dimensions, "
      << "and generated " << vegas_state_->bins_max << " bins.\n\t"
      << "Integration volume: " << vegas_state_->vol << ".";
    grid_->r_boxes = std::pow( vegas_state_->bins, function_->dim );

    result_ = result;
    err_result_ = abserr;
  }

  void
  IntegratorVegas::warmup( std::vector<double>& x_low, std::vector<double>& x_up, unsigned int ncall )
  {
    // start by preparing the grid/state
    vegas_state_.reset( gsl_monte_vegas_alloc( function_->dim ) );
    gsl_monte_vegas_params_set( vegas_state_.get(), &vegas_params_ );
    // then perform a first integration with the given calls count
    double result = 0., abserr = 0.;
    const int res = gsl_monte_vegas_integrate( function_.get(),
      &x_low[0], &x_up[0],
      function_->dim, ncall, rng_.get(), vegas_state_.get(),
      &result, &abserr );
    // ensure the operation was successful
    if ( res != GSL_SUCCESS )
      throw CG_ERROR( "Integrator:vegas" )
        << "Failed to warm-up the Vegas grid.\n\t"
        << "GSL error: " << gsl_strerror( res ) << ".";
    CG_INFO( "Integrator:vegas" )
      << "Finished the Vegas warm-up.";
  }

  double
  IntegratorVegas::eval( const std::vector<double>& x )
  {
    if ( !input_params_->generation().treat )
      return function_->f( (double*)&x[0], function_->dim, (void*)&input_params_ );
    //--- treatment of the integration grid
    double w = grid_->r_boxes;
    std::vector<double> x_new( x.size() );
    for ( unsigned short j = 0; j < function_->dim; ++j ) {
      //--- find surrounding coordinates and interpolate
      const double z = x[j]*vegas_state_->bins;
      const unsigned int id = z; // coordinate of point before
      const double rel_pos = z-id; // position between coordinates (norm.)
      const double bin_width = ( id == 0 )
        ? COORD( vegas_state_, 1, j )
        : COORD( vegas_state_, id+1, j )-COORD( vegas_state_, id, j );
      //--- build new coordinate from linear interpolation
      x_new[j] = COORD( vegas_state_, id+1, j )-bin_width*( 1.-rel_pos );
      w *= bin_width;
    }
    return w*function_->f( (double*)&x_new[0], function_->dim, (void*)&input_params_ );
  }

  std::ostream&
  operator<<( std::ostream& os, const IntegratorVegas::Mode& mode )
  {
    switch ( mode ) {
      case IntegratorVegas::Mode::importance:
        return os << "importance";
      case IntegratorVegas::Mode::importanceOnly:
        return os << "importance-only";
      case IntegratorVegas::Mode::stratified:
        return os << "stratified";
    }
    return os;
  }
}

REGISTER_INTEGRATOR( "Vegas", IntegratorVegas )
