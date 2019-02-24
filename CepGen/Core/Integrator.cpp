#include "CepGen/Core/Integrator.h"
#include "CepGen/Core/GridParameters.h"
#include "CepGen/Core/utils.h"
#include "CepGen/Core/Exception.h"

#include "CepGen/Parameters.h"
#include "CepGen/Processes/GenericProcess.h"
#include "CepGen/Hadronisers/GenericHadroniser.h"
#include "CepGen/Event/Event.h"

#include <thread>
#include <math.h>

#include <gsl/gsl_monte_miser.h>
#define COORD(s,i,j) ((s)->xi[(i)*(s)->dim + (j)])

namespace cepgen
{
  Integrator::Integrator( unsigned int ndim, double integrand( double*, size_t, void* ), Parameters* params ) :
    ps_bin_( INVALID_BIN ), input_params_( params, []( Parameters* ){} ),
    function_( new gsl_monte_function{ integrand, ndim, (void*)input_params_.get() } ),
    rng_( gsl_rng_alloc( input_params_->integration().rng_engine ), gsl_rng_free ),
    grid_( new GridParameters( ndim ) )
  {
    //--- initialise the random number generator

    unsigned long seed = ( input_params_->integration().rng_seed > 0 )
      ? input_params_->integration().rng_seed
      : time( nullptr ); // seed with time
    gsl_rng_set( rng_.get(), seed );

    //--- a bit of printout for debugging

    CG_DEBUG( "Integrator:build" )
      << "Number of integration dimensions: " << function_->dim << ",\n\t"
      << "Number of function calls:         " << input_params_->integration().ncvg << ",\n\t"
      << "Random numbers generator:         " << gsl_rng_name( rng_.get() ) << ".";
    switch ( input_params_->integration().type ) {
      case IntegratorType::Vegas:
        CG_DEBUG( "Integrator:build" ) << "Vegas parameters:\n\t"
          << "Number of iterations in Vegas: " << input_params_->integration().vegas.iterations << ",\n\t"
          << "α-value: " << input_params_->integration().vegas.alpha << ",\n\t"
          << "Verbosity: " << input_params_->integration().vegas.verbose << ",\n\t"
          << "Grid interpolation mode: " << (Integrator::VegasMode)input_params_->integration().vegas.mode << ".";
        break;
      case IntegratorType::MISER:
        CG_DEBUG( "Integrator:build" ) << "MISER parameters:\n\t"
          << "Number of calls: " << input_params_->integration().miser.min_calls << ", "
          << "per bisection: " << input_params_->integration().miser.min_calls_per_bisection << ",\n\t"
          << "Estimate fraction: " << input_params_->integration().miser.estimate_frac << ",\n\t"
          << "α-value: " << input_params_->integration().miser.alpha << ",\n\t"
          << "Dither: " << input_params_->integration().miser.dither << ".";
        break;
      case IntegratorType::plain:
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

    //--- integration bounds
    std::vector<double> x_low( function_->dim, 0. ), x_up( function_->dim, 1. );

    //--- launch integration
    switch ( input_params_->integration().type ) {
      case IntegratorType::plain: {
        std::unique_ptr<gsl_monte_plain_state,void(*)( gsl_monte_plain_state* )>
          pln_state( gsl_monte_plain_alloc( function_->dim ), gsl_monte_plain_free );
        res = gsl_monte_plain_integrate( function_.get(),
          &x_low[0], &x_up[0],
          function_->dim, input_params_->integration().ncvg,
          rng_.get(), pln_state.get(),
          &result, &abserr );
      } break;
      case IntegratorType::Vegas: {
        //----- warmup (prepare the grid)
        res = warmupVegas( x_low, x_up, 25000 );
        //----- integration
        unsigned short it_chisq = 0;
        do {
          res = gsl_monte_vegas_integrate( function_.get(),
            &x_low[0], &x_up[0],
            function_->dim, 0.2 * input_params_->integration().ncvg,
            rng_.get(), veg_state_.get(),
            &result, &abserr );
          CG_LOG( "Integrator:integrate" )
            << "\t>> at call " << ( ++it_chisq ) << ": "
            << Form( "average = %10.6f   "
                     "sigma = %10.6f   chi2 = %4.3f.",
                     result, abserr,
                     gsl_monte_vegas_chisq( veg_state_.get() ) );
        } while ( fabs( gsl_monte_vegas_chisq( veg_state_.get() )-1. )
                > input_params_->integration().vegas_chisq_cut-1. );
        CG_DEBUG( "Integrator:integrate" )
          << "Vegas grid information:\n\t"
          << "ran for " << veg_state_->dim << " dimensions, and generated " << veg_state_->bins_max << " bins.\n\t"
          << "Integration volume: " << veg_state_->vol << ".";
        grid_->r_boxes = std::pow( veg_state_->bins, function_->dim );
      } break;
      case IntegratorType::MISER: {
        std::unique_ptr<gsl_monte_miser_state,void(*)( gsl_monte_miser_state* )>
          mis_state( gsl_monte_miser_alloc( function_->dim ), gsl_monte_miser_free );
        gsl_monte_miser_params_set( mis_state.get(), &input_params_->integration().miser );
        res = gsl_monte_miser_integrate( function_.get(),
          &x_low[0], &x_up[0],
          function_->dim, input_params_->integration().ncvg,
          rng_.get(), mis_state.get(),
          &result, &abserr );
      } break;
    }

    input_params_->integration().result = result;
    input_params_->integration().err_result = abserr;

    if ( input_params_->hadroniser() )
      input_params_->hadroniser()->setCrossSection( result, abserr );

    return res;
  }

  int
  Integrator::warmupVegas( std::vector<double>& x_low, std::vector<double>& x_up, unsigned int ncall )
  {
    // start by preparing the grid/state
    veg_state_.reset( gsl_monte_vegas_alloc( function_->dim ) );
    gsl_monte_vegas_params_set( veg_state_.get(), &input_params_->integration().vegas );
    // then perform a first integration with the given calls count
    double result = 0., abserr = 0.;
    int res = gsl_monte_vegas_integrate( function_.get(),
      &x_low[0], &x_up[0],
      function_->dim, ncall,
      rng_.get(), veg_state_.get(),
      &result, &abserr );
    // ensure the operation was successful
    if ( res != GSL_SUCCESS )
      throw CG_ERROR( "Integrator:vegas" )
        << "Failed to warm-up the Vegas grid.\n\t"
        << "GSL error: " << gsl_strerror( res ) << ".";
    CG_INFO( "Integrator:vegas" )
      << "Finished the Vegas warm-up.";
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

    if ( ps_bin_ != INVALID_BIN ) {
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
        x[i] = ( uniform() + grid_n[i] ) * GridParameters::INV_M_BIN;
      // get weight for selected x value
      weight = eval( x );
      if ( weight <= 0. )
        continue;
      if ( weight > y )
        break;
    }

    if ( weight <= grid_->f_max[ps_bin_] )
      ps_bin_ = INVALID_BIN;
    else {
      //--- if weight is higher than local or global maximum,
      //    init correction cycle
      grid_->f_max_old = grid_->f_max[ps_bin_];
      grid_->f_max[ps_bin_] = weight;
      grid_->f_max_diff = weight-grid_->f_max_old;
      if ( weight > grid_->f_max_global )
        grid_->f_max_global = weight;
      grid_->correc = ( grid_->num[ps_bin_]-1. ) * grid_->f_max_diff / grid_->f_max_global - 1.;

      CG_DEBUG("Integrator::generateOne")
        << "Correction " << grid_->correc << " will be applied for phase space bin " << ps_bin_ << ".";
    }

    // return with an accepted event
    if ( weight > 0. )
      storeEvent( x, callback );
  }

  void
  Integrator::generate( unsigned long num_events, std::function<void( const Event&, unsigned long )> callback )
  {
    if ( num_events < 1 )
      num_events = input_params_->generation().maxgen;
    try {
      while ( input_params_->numGeneratedEvents() < num_events )
        generateOne( callback );
    } catch ( const Exception& e ) { return; }
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
        xtmp[k] = ( uniform() + grid_n[k] ) * GridParameters::INV_M_BIN;
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
    //--- start by computing the matrix element for that point
    const double weight = eval( x );

    //--- reject if unphysical
    if ( weight <= 0. )
      return false;

    {
      if ( input_params_->numGeneratedEvents() % input_params_->generation().gen_print_every == 0 ) {
        CG_INFO( "Integrator:store" )
          << "Generated events: " << input_params_->numGeneratedEvents();
        input_params_->process()->last_event->dump();
      }
      const Event& last_event = *input_params_->process()->last_event;
      if ( callback )
        callback( last_event, input_params_->numGeneratedEvents() );
      input_params_->addGenerationTime( last_event.time_total );
    }
    return true;
  }

  //-----------------------------------------------------------------------------------------------
  // initial preparation run before the generation of unweighted events
  //-----------------------------------------------------------------------------------------------

  void
  Integrator::computeGenerationParameters()
  {
    input_params_->setStorage( false );

    if ( input_params_->generation().treat
      && input_params_->integration().type != IntegratorType::Vegas ) {
      CG_INFO( "Integrator:setGen" )
        << "Treat switched on without a proper Vegas grid; running a warm-up beforehand.";
      std::vector<double> x_low( function_->dim, 0. ), x_up( function_->dim, 1. );
      int res = warmupVegas( x_low, x_up, 25000 );
      if ( res != GSL_SUCCESS )
        throw CG_FATAL( "Integrator::setGen" )
          << "Failed to perform a Vegas warm-up.\n\t"
          << "Try to re-run while disabling integrand treatment...";
    }
    CG_INFO( "Integrator:setGen" )
      << "Preparing the grid (" << input_params_->generation().num_points << " points/bin) "
      << "for the generation of unweighted events.";

    const double inv_num_points = 1./input_params_->generation().num_points;
    std::vector<double> x( function_->dim, 0. );
    std::vector<unsigned short> n( function_->dim, 0 );;

    // ...
    double sum = 0., sum2 = 0., sum2p = 0.;

    //--- main loop
    for ( unsigned int i = 0; i < grid_->max; ++i ) {
      unsigned int jj = i;
      for ( unsigned int j = 0; j < function_->dim; ++j ) {
        unsigned int tmp = jj*GridParameters::INV_M_BIN;
        n[j] = jj-tmp*GridParameters::M_BIN;
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
      for ( unsigned int j = 0; j < input_params_->generation().num_points; ++j ) {
        for ( unsigned int k = 0; k < function_->dim; ++k )
          x[k] = ( uniform()+n[k] ) * GridParameters::INV_M_BIN;
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
    input_params_->setStorage( true );
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
    if ( !input_params_->generation().treat )
      return function_->f( (double*)&x[0], function_->dim, (void*)input_params_.get() );
    //--- treatment of the integration grid
    double w = grid_->r_boxes;
    std::vector<double> x_new( x.size() );
    for ( unsigned short j = 0; j < function_->dim; ++j ) {
      //--- find surrounding coordinates and interpolate
      const double z = x[j]*veg_state_->bins;
      const unsigned int id = z; // coordinate of point before
      const double rel_pos = z-id; // position between coordinates (norm.)
      const double bin_width = ( id == 0 )
        ? COORD( veg_state_, 1, j )
        : COORD( veg_state_, id+1, j )-COORD( veg_state_, id, j );
      //--- build new coordinate from linear interpolation
      x_new[j] = COORD( veg_state_, id+1, j )-bin_width*( 1.-rel_pos );
      w *= bin_width;
    }
    return w*function_->f( (double*)&x_new[0], function_->dim, (void*)input_params_.get() );
  }

  double
  Integrator::uniform() const
  {
    return gsl_rng_uniform( rng_.get() );
  }

  //------------------------------------------------------------------------------------------------

  std::ostream&
  operator<<( std::ostream& os, const IntegratorType& type )
  {
    switch ( type ) {
      case IntegratorType::plain:
        return os << "plain";
      case IntegratorType::Vegas:
        return os << "Vegas";
      case IntegratorType::MISER:
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

