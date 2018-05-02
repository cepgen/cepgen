#include "CepGen/Core/Integrator.h"
#include "CepGen/Core/GridParameters.h"
#include "CepGen/Core/utils.h"
#include "CepGen/Core/Exception.h"

#include "CepGen/Parameters.h"
#include "CepGen/Processes/GenericProcess.h"
#include "CepGen/Hadronisers/GenericHadroniser.h"

#include <thread>
#include <math.h>

#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_monte_miser.h>

namespace CepGen
{
  Integrator::Integrator( unsigned int ndim, double integrand( double*, size_t, void* ), Parameters* params ) :
    ps_bin_( 0 ), input_params_( params ),
    function_( new gsl_monte_function{ integrand, ndim, (void*)input_params_ } ),
    rng_( gsl_rng_alloc( input_params_->integrator.rng_engine ), gsl_rng_free ),
    grid_( new GridParameters )
  {
    //--- initialise the random number generator

    unsigned long seed = ( input_params_->integrator.rng_seed > 0 )
      ? input_params_->integrator.rng_seed
      : time( nullptr ); // seed with time
    gsl_rng_set( rng_.get(), seed );

    //--- a bit of printout for debugging

    CG_DEBUG( "Integrator:build" )
      << "Number of integration dimensions: " << function_->dim << ",\n\t"
      << "Number of function calls:         " << input_params_->integrator.ncvg << ",\n\t"
      << "Random numbers generator:         " << gsl_rng_name( rng_.get() ) << ".";
    switch ( input_params_->integrator.type ) {
      case Type::Vegas:
        CG_DEBUG( "Integrator:build" ) << "Vegas parameters:\n\t"
          << "Number of iterations in Vegas: " << input_params_->integrator.vegas.iterations << ",\n\t"
          << "α-value: " << input_params_->integrator.vegas.alpha << ",\n\t"
          << "Verbosity: " << input_params_->integrator.vegas.verbose << ",\n\t"
          << "Grid interpolation mode: " << (Integrator::VegasMode)input_params_->integrator.vegas.mode << ".";
        break;
      case Type::MISER:
        CG_DEBUG( "Integrator:build" ) << "MISER parameters:\n\t"
          << "Number of calls: " << input_params_->integrator.miser.min_calls << ", "
          << "per bisection: " << input_params_->integrator.miser.min_calls_per_bisection << ",\n\t"
          << "Estimate fraction: " << input_params_->integrator.miser.estimate_frac << ",\n\t"
          << "α-value: " << input_params_->integrator.miser.alpha << ",\n\t"
          << "Dither: " << input_params_->integrator.miser.dither << ".";
        break;
      case Type::plain:
        break;
    }
  }

  Integrator::~Integrator()
  {}

  //-----------------------------------------------------------------------------------------------
  // cross section computation part
  //-----------------------------------------------------------------------------------------------

  int
  Integrator::integrate( double& result, double& abserr )
  {
    int res = -1;
    gsl_monte_plain_state* pln_state = nullptr;
    gsl_monte_vegas_state* veg_state = nullptr;
    gsl_monte_miser_state* mis_state = nullptr;
    const Integrator::Type algorithm = input_params_->integrator.type;

    //--- integration bounds
    std::vector<double> x_low( function_->dim, 0. ), x_up( function_->dim, 1. );

    //--- prepare integrator
    if ( algorithm == Type::plain )
      pln_state = gsl_monte_plain_alloc( function_->dim );
    else if ( algorithm == Type::Vegas )
      veg_state = gsl_monte_vegas_alloc( function_->dim );
    else if ( algorithm == Type::MISER )
      mis_state = gsl_monte_miser_alloc( function_->dim );

    if ( algorithm == Type::plain )
      res = gsl_monte_plain_integrate( function_.get(),
        &x_low[0], &x_up[0],
        function_->dim, input_params_->integrator.ncvg,
        rng_.get(), pln_state,
        &result, &abserr );
    else if ( algorithm == Type::Vegas ) {
      gsl_monte_vegas_params_set( veg_state, &input_params_->integrator.vegas );
      //----- Vegas warmup (prepare the grid)
      res = gsl_monte_vegas_integrate( function_.get(),
        &x_low[0], &x_up[0],
        function_->dim, 25000,
        rng_.get(), veg_state,
        &result, &abserr );
      CG_INFO( "Integrator:integrate" )
        << "Finished the Vegas warm-up.\n\t"
        << "Will now iterate until χ² < " << input_params_->integrator.vegas_chisq_cut << ".";
      //----- integration
      unsigned short it_chisq = 0;
      do {
        res = gsl_monte_vegas_integrate( function_.get(),
          &x_low[0], &x_up[0],
          function_->dim, 0.2 * input_params_->integrator.ncvg,
          rng_.get(), veg_state,
          &result, &abserr );
        CG_LOG( "Integrator:integrate" )
          << "\t>> at call " << ( it_chisq+1 ) << ": "
          << Form( "average = %10.6f   "
                   "sigma = %10.6f   chi2 = %4.3f.",
                   result, abserr,
                   gsl_monte_vegas_chisq( veg_state ) );
        it_chisq++;
      } while ( fabs( gsl_monte_vegas_chisq( veg_state )-1. )
              > input_params_->integrator.vegas_chisq_cut-1. );
    }
    //----- integration
    else if ( algorithm == Type::MISER ) {
      gsl_monte_miser_params_set( mis_state, &input_params_->integrator.miser );
      res = gsl_monte_miser_integrate( function_.get(),
        &x_low[0], &x_up[0],
        function_->dim, input_params_->integrator.ncvg,
        rng_.get(), mis_state,
        &result, &abserr );
    }

    //--- clean integrator
    if ( algorithm == Type::plain )
      gsl_monte_plain_free( pln_state );
    else if ( algorithm == Type::Vegas ) {
      CG_DEBUG( "Integrator:integrate" )
        << "Vegas grid information:\n\t"
        << "ran for " << veg_state->dim << " dimensions, and generated " << veg_state->bins_max << " bins.\n\t"
        << "Integration volume: " << veg_state->vol << ".";
      gsl_monte_vegas_free( veg_state );
    }
    else if ( algorithm == Type::MISER )
      gsl_monte_miser_free( mis_state );

    input_params_->integrator.result = result;
    input_params_->integrator.err_result = abserr;

    if ( input_params_->hadroniser() )
      input_params_->hadroniser()->setCrossSection( result, abserr );

    return res;
  }

  //-----------------------------------------------------------------------------------------------
  // events generation part
  //-----------------------------------------------------------------------------------------------

  void
  Integrator::generateOne( std::function<void( const Event&, unsigned long )> callback )
  {
    if ( !grid_->gen_prepared )
      computeGenerationParameters();

    std::vector<double> x( function_->dim, 0. );

    //--- correction cycles

    if ( ps_bin_ != 0 ) {
      bool has_correction = false;
      while ( !correctionCycle( x, has_correction ) ) {}
      if ( has_correction ) {
        storeEvent( x, callback );
        return;
      }
    }

    double weight = 0.;

    //--- normal generation cycle

    while ( true ) {
      double y = -1.;
      //----- select a and reject if fmax is too small
      while ( true ) {
        // ...
        ps_bin_ = uniform() * grid_->max;
        y = uniform() * grid_->f_max_global;
        grid_->num[ps_bin_] += 1;
        if ( y <= grid_->f_max[ps_bin_] )
          break;
      }
      // shoot a point x in this bin
      std::vector<unsigned short> grid_n = grid_->n_map.at( ps_bin_ );
      for ( unsigned int i = 0; i < function_->dim; ++i )
        x[i] = ( uniform() + grid_n[i] ) * GridParameters::inv_mbin_;
      // get weight for selected x value
      weight = eval( x );
      if ( weight <= 0. )
        continue;
      if ( weight > y )
        break;
    }

    if ( weight <= grid_->f_max[ps_bin_] )
      ps_bin_ = 0;
    // init correction cycle if weight is higher than fmax or ffmax
    else if ( weight <= grid_->f_max_global ) {
      grid_->f_max_old = grid_->f_max[ps_bin_];
      grid_->f_max[ps_bin_] = weight;
      grid_->f_max_diff = weight-grid_->f_max_old;
      grid_->correc = ( grid_->num[ps_bin_]-1. ) * grid_->f_max_diff / grid_->f_max_global - 1.;
    }
    else {
      grid_->f_max_old = grid_->f_max[ps_bin_];
      grid_->f_max[ps_bin_] = weight;
      grid_->f_max_diff = weight-grid_->f_max_old;
      grid_->f_max_global = weight;
      grid_->correc = ( grid_->num[ps_bin_] - 1. ) * grid_->f_max_diff / grid_->f_max_global * weight / grid_->f_max_global - 1.;
    }

    CG_DEBUG("Integrator::generateOne")
      << "Correction " << grid_->correc << " will be applied for phase space bin " << ps_bin_ << ".";

    // return with an accepted event
    if ( weight > 0. )
      storeEvent( x, callback );
  }

  void
  Integrator::generate( unsigned long num_events, std::function<void( const Event&, unsigned long )> callback, const Timer* tmr )
  {
    if ( num_events < 1 )
      num_events = input_params_->generation.maxgen;
    while ( input_params_->generation.ngen < num_events )
      generateOne( callback );
  }

  bool
  Integrator::correctionCycle( std::vector<double>& x, bool& has_correction )
  {
    CG_DEBUG_LOOP( "Integrator:correction" )
      << "Correction cycles are started.\n\t"
      << "bin = " << ps_bin_ << "\t"
      << "correc = " << grid_->correc << "\t"
      << "corre2 = " << grid_->correc2 << ".";

    if ( grid_->correc >= 1. )
      grid_->correc -= 1.;

    if ( uniform() < grid_->correc ) {
      grid_->correc = -1.;
      std::vector<double> xtmp( function_->dim );
      // Select x values in phase space bin
      const std::vector<unsigned short> grid_n = grid_->n_map.at( ps_bin_ );
      for ( unsigned int k = 0; k < function_->dim; ++k )
        xtmp[k] = ( uniform() + grid_n[k] ) * GridParameters::inv_mbin_;
      const double weight = eval( xtmp );
      // Parameter for correction of correction
      if ( weight > grid_->f_max[ps_bin_] ) {
        grid_->f_max2 = std::max( grid_->f_max2, weight );
        grid_->correc += 1.;
        grid_->correc2 -= 1.;
      }
      // Accept event
      if ( weight >= grid_->f_max_diff*uniform() + grid_->f_max_old ) {
        x = xtmp;
        has_correction = true;
        return true;
      }
      return false;
    }
    // Correction if too big weight is found while correction
    // (All your bases are belong to us...)
    if ( grid_->f_max2 > grid_->f_max[ps_bin_] ) {
      grid_->f_max_old = grid_->f_max[ps_bin_];
      grid_->f_max[ps_bin_] = grid_->f_max2;
      grid_->f_max_diff = grid_->f_max2-grid_->f_max_old;
      grid_->correc = ( grid_->num[ps_bin_]-1. ) * grid_->f_max_diff / grid_->f_max_global;
      if ( grid_->f_max2 >= grid_->f_max_global ) {
        grid_->correc *= grid_->f_max2 / grid_->f_max_global;
        grid_->f_max_global = grid_->f_max2;
      }
      grid_->correc -= grid_->correc2;
      grid_->correc2 = 0.;
      grid_->f_max2 = 0.;
      return false;
    }
    return true;
  }

  bool
  Integrator::storeEvent( const std::vector<double>& x, std::function<void( const Event&, unsigned long )> callback )
  {
    input_params_->setStorage( true );
    const double weight = eval( x );
    input_params_->setStorage( false );

    if ( weight <= 0. )
      return false;

    {
      if ( input_params_->generation.ngen % input_params_->generation.gen_print_every == 0 ) {
        CG_INFO( "Integrator:store" )
          << "Generated events: " << input_params_->generation.ngen;
        input_params_->process()->last_event->dump();
      }
      input_params_->generation.ngen += 1;

      if ( callback )
        callback( *input_params_->process()->last_event, input_params_->generation.ngen );
    }
    return true;
  }

  //-----------------------------------------------------------------------------------------------
  // initial preparation run before the generation of unweighted events
  //-----------------------------------------------------------------------------------------------

  void
  Integrator::computeGenerationParameters()
  {
    input_params_->generation.ngen = 0;
    input_params_->setStorage( false );

    CG_INFO( "Integrator:setGen" )
      << "Preparing the grid (" << input_params_->generation.num_points << " points/bin) "
      << "for the generation of unweighted events.";

    grid_->max = pow( GridParameters::mbin_, function_->dim );
    const double inv_num_points = 1./input_params_->generation.num_points;

    if ( function_->dim > GridParameters::max_dimensions_ )
      throw CG_FATAL( "Integrator:setGen" )
        << "Number of dimensions to integrate exceeds the maximum number, "
        << GridParameters::max_dimensions_ << ".";

    grid_->f_max = std::vector<double>( grid_->max, 0. );
    grid_->num.reserve( grid_->max );

    std::vector<double> x( function_->dim, 0. );
    std::vector<unsigned short> n( function_->dim, 0 );;

    // ...
    double sum = 0., sum2 = 0., sum2p = 0.;

    //--- main loop
    for ( unsigned int i = 0; i < grid_->max; ++i ) {
      unsigned int jj = i;
      for ( unsigned int j = 0; j < function_->dim; ++j ) {
        unsigned int tmp = jj*GridParameters::inv_mbin_;
        n[j] = jj-tmp*GridParameters::mbin_;
        jj = tmp;
      }
      grid_->n_map[i] = n;
      if ( CG_EXCEPT_MATCH( "Integrator:setGen", debugInsideLoop ) ) {
        std::ostringstream os;
        for ( const auto& ni : n )
          os << ni << " ";
        CG_DEBUG_LOOP( "Integrator:setGen" )
          << "n-vector for bin " << i << ": " << os.str();
      }
      double fsum = 0., fsum2 = 0.;
      for ( unsigned int j = 0; j < input_params_->generation.num_points; ++j ) {
        for ( unsigned int k = 0; k < function_->dim; ++k )
          x[k] = ( uniform()+n[k] ) * GridParameters::inv_mbin_;
        const double weight = eval( x );
        grid_->f_max[i] = std::max( grid_->f_max[i], weight );
        fsum += weight;
        fsum2 += weight*weight;
      }
      const double av = fsum*inv_num_points, av2 = fsum2*inv_num_points, sig2 = av2-av*av;
      sum += av;
      sum2 += av2;
      sum2p += sig2;
      grid_->f_max_global = std::max( grid_->f_max_global, grid_->f_max[i] );

      // per-bin debugging loop
      if ( CG_EXCEPT_MATCH( "Integrator:setGen", debugInsideLoop ) ) {
        const double sig = sqrt( sig2 );
        const double eff = ( grid_->f_max[i] != 0. )
          ? grid_->f_max[i]/av
          : 1.e4;
        std::ostringstream os;
        for ( unsigned int j = 0; j < function_->dim; ++j )
          os << ( j != 0 ? ", " : "" ) << n[j];
        CG_DEBUG_LOOP( "Integrator:setGen" )
          << "In iteration #" << i << ":\n\t"
          << "av   = " << av << "\n\t"
          << "sig  = " << sig << "\n\t"
          << "fmax = " << grid_->f_max[i] << "\n\t"
          << "eff  = " << eff << "\n\t"
          << "n = (" << os.str() << ")";
      }
    } // end of main loop

    const double inv_max = 1./grid_->max;
    sum *= inv_max;
    sum2 *= inv_max;
    sum2p *= inv_max;

    const double sig = sqrt( sum2-sum*sum ), sigp = sqrt( sum2p );

    double eff1 = 0.;
    for ( unsigned int i = 0; i < grid_->max; ++i )
      eff1 += sum/grid_->max*grid_->f_max[i];
    const double eff2 = sum/grid_->f_max_global;

    CG_DEBUG( "Integrator:setGen" )
      << "Average function value         = sum   = " << sum << "\n\t"
      << "Average squared function value = sum2  = " << sum2 << "\n\t"
      << "Overall standard deviation     = sig   = " << sig << "\n\t"
      << "Average standard deviation     = sigp  = " << sigp << "\n\t"
      << "Maximum function value         = f_max = " << grid_->f_max_global << "\n\t"
      << "Average inefficiency           = eff1  = " << eff1 << "\n\t"
      << "Overall inefficiency           = eff2  = " << eff2;

    grid_->gen_prepared = true;
    CG_INFO( "Integrator:setGen" ) << "Grid prepared! Now launching the production.";
  }

  //------------------------------------------------------------------------------------------------
  // helper / alias methods
  //------------------------------------------------------------------------------------------------

  unsigned short
  Integrator::dimensions() const
  {
    if ( !function_ )
      return 0;
    return function_->dim;
  }

  double
  Integrator::eval( const std::vector<double>& x )
  {
    return function_->f( (double*)&x[0], function_->dim, (void*)input_params_ );
  }

  double
  Integrator::uniform() const
  {
    return gsl_rng_uniform( rng_.get() );
  }

  //------------------------------------------------------------------------------------------------

  std::ostream&
  operator<<( std::ostream& os, const Integrator::Type& type )
  {
    switch ( type ) {
      case Integrator::Type::plain:
        return os << "plain";
      case Integrator::Type::Vegas:
        return os << "Vegas";
      case Integrator::Type::MISER:
        return os << "MISER";
    }
    return os;
  }

  std::ostream&
  operator<<( std::ostream& os, const Integrator::VegasMode& mode )
  {
    switch ( mode ) {
      case Integrator::VegasMode::importance:
        return os << "importance";
      case Integrator::VegasMode::importanceOnly:
        return os << "importance-only";
      case Integrator::VegasMode::stratified:
        return os << "stratified";
    }
    return os;
  }
}

