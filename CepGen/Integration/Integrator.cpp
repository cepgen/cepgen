#include "CepGen/Integration/Integrator.h"
#include "CepGen/Integration/GridParameters.h"

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
    name_( params.name<std::string>() ),
    ncvg_( params.get<int>( "numFunctionCalls", 50000 ) ),
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
    grid_.reset( new GridParameters( ndim ) );
    input_params_ = &params;
    function_.reset( new gsl_monte_function{ integrand, ndim, (void*)input_params_ } );

    CG_DEBUG( "Integrator:function" )
      << "Number of integration dimensions: " << function_->dim << ".";

    initialised_ = false;
  }

  Integrator::~Integrator()
  {}

  //-----------------------------------------------------------------------------------------------
  // events generation part
  //-----------------------------------------------------------------------------------------------

  void
  Integrator::generateOne( Event::callback callback )
  {
    if ( !grid_->gen_prepared )
      computeGenerationParameters();

    std::vector<double> xtmp( function_->dim );

    //--- correction cycles

    if ( ps_bin_ != INVALID_BIN ) {
      bool has_correction = false;
      while ( !correctionCycle( xtmp, has_correction ) ) {}
      if ( has_correction ) {
        storeEvent( xtmp, callback );
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
        ps_bin_ = uniform() * grid_->size();
        y = uniform() * grid_->globalMax();
        grid_->setTrial( ps_bin_ );
        if ( y <= grid_->maxValue( ps_bin_ ) )
          break;
      }
      // shoot a point x in this bin
      grid_->shoot( rng_.get(), ps_bin_, xtmp );
      // get weight for selected x value
      weight = eval( xtmp );
      if ( weight <= 0. )
        continue;
      if ( weight > y )
        break;
    }

    if ( weight <= grid_->maxValue( ps_bin_ ) )
      ps_bin_ = INVALID_BIN;
    else {
      //--- if weight is higher than local or global maximum,
      //    init correction cycle
      grid_->f_max_old = grid_->maxValue( ps_bin_ );
      grid_->f_max_diff = weight-grid_->f_max_old;
      grid_->setValue( ps_bin_, weight );
      grid_->correc = ( grid_->numPoints( ps_bin_ )-1. ) * grid_->f_max_diff / grid_->globalMax() - 1.;

      CG_DEBUG("Integrator:generateOne")
        << "Correction " << grid_->correc << " will be applied for phase space bin " << ps_bin_ << ".";
    }

    // return with an accepted event
    if ( weight > 0. )
      storeEvent( xtmp, callback );
  }

  void
  Integrator::generate( unsigned long num_events, Event::callback callback )
  {
    if ( num_events < 1 )
      num_events = input_params_->generation().maxgen;
    for ( auto& mod : input_params_->outputModulesSequence() )
      mod->initialise( *input_params_ );
    try {
      while ( input_params_->numGeneratedEvents() < num_events )
        generateOne( callback );
    } catch ( const Exception& ) { return; }
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
      grid_->shoot( rng_.get(), ps_bin_, xtmp );
      const double weight = eval( xtmp );
      // Parameter for correction of correction
      if ( weight > grid_->maxValue( ps_bin_ ) ) {
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
    if ( grid_->f_max2 > grid_->maxValue( ps_bin_ ) ) {
      grid_->f_max_old = grid_->maxValue( ps_bin_ );
      grid_->f_max_diff = grid_->f_max2-grid_->f_max_old;
      grid_->correc = ( grid_->numPoints( ps_bin_ )-1. ) * grid_->f_max_diff / grid_->globalMax();
      if ( grid_->f_max2 >= grid_->globalMax() )
        grid_->correc *= grid_->f_max2 / grid_->globalMax();
      grid_->setValue( ps_bin_, grid_->f_max2 );
      grid_->correc -= grid_->correc2;
      grid_->correc2 = 0.;
      grid_->f_max2 = 0.;
      return false;
    }
    return true;
  }

  bool
  Integrator::storeEvent( const std::vector<double>& x, Event::callback callback )
  {
    //--- start by computing the matrix element for that point
    const double weight = eval( x );

    //--- reject if unphysical
    if ( weight <= 0. )
      return false;

    const auto ngen = input_params_->numGeneratedEvents();
    if ( input_params_->process().hasEvent() ) {
      auto& event = input_params_->process().event();
      if ( ( ngen+1 ) % input_params_->generation().gen_print_every == 0 )
        CG_INFO( "Integrator:store" )
          << "Generated events: " << ngen+1 << "\n"
          << event;
      if ( callback )
        callback( event, ngen );
      for ( auto& mod : input_params_->outputModulesSequence() )
        *mod << event;
      input_params_->addGenerationTime( event.time_total );
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

    if ( input_params_->generation().treat && name_ != "Vegas" ) {
      throw CG_FATAL( "Integrator:setGen" )
        << "Treat switched on without a proper Vegas grid!\n\t"
        << "Try to re-run while disabling integrand treatment...";
    }
    CG_INFO( "Integrator:setGen" )
      << "Preparing the grid (" << input_params_->generation().num_points << " points/bin) "
      << "for the generation of unweighted events.";

    const double inv_num_points = 1./input_params_->generation().num_points;
    std::vector<double> point_coord( function_->dim, 0. );
    if ( point_coord.size() < grid_->n( 0 ).size() )
      throw CG_FATAL( "GridParameters:shoot" )
        << "Coordinates vector multiplicity is insufficient!";

    // ...
    double sum = 0., sum2 = 0., sum2p = 0.;

    utils::ProgressBar prog_bar( grid_->size(), 5 );

    //--- main loop
    for ( unsigned int i = 0; i < grid_->size(); ++i ) {
      double fsum = 0., fsum2 = 0.;
      for ( unsigned int j = 0; j < input_params_->generation().num_points; ++j ) {
        grid_->shoot( rng_.get(), i, point_coord );
        const double weight = eval( point_coord );
        grid_->setValue( i, weight );
        fsum += weight;
        fsum2 += weight*weight;
      }
      const double av = fsum*inv_num_points, av2 = fsum2*inv_num_points;
      const double sig2 = av2-av*av;
      sum += av;
      sum2 += av2;
      sum2p += sig2;

      // per-bin debugging loop
      if ( CG_LOG_MATCH( "Integrator:setGen", debugInsideLoop ) ) {
        const double sig = sqrt( sig2 );
        const double eff = ( grid_->maxValue( i ) != 0. )
          ? av/grid_->maxValue( i )
          : 0.;
        CG_DEBUG_LOOP( "Integrator:setGen" )
          << "n-vector for bin " << i << ": " << utils::repr( grid_->n( i ) ) << "\n\t"
          << "av   = " << av << "\n\t"
          << "sig  = " << sig << "\n\t"
          << "fmax = " << grid_->maxValue( i ) << "\n\t"
          << "eff  = " << eff;
      }
      prog_bar.update( i+1 );
    } // end of main loop

    const double inv_max = 1./grid_->size();
    sum *= inv_max;
    sum2 *= inv_max;
    sum2p *= inv_max;

    const double sig = sqrt( sum2-sum*sum ), sigp = sqrt( sum2p );

    double eff1 = 0.;
    for ( unsigned int i = 0; i < grid_->size(); ++i )
      eff1 += sum/grid_->size()*grid_->maxValue( i );
    const double eff2 = sum/grid_->globalMax();

    CG_DEBUG( "Integrator:setGen" )
      << "Average function value         = " << sum << "\n\t"
      << "Average squared function value = " << sum2 << "\n\t"
      << "Overall standard deviation     = " << sig << "\n\t"
      << "Average standard deviation     = " << sigp << "\n\t"
      << "Maximum function value         = " << grid_->globalMax() << "\n\t"
      << "Average inefficiency           = " << eff1 << "\n\t"
      << "Overall inefficiency           = " << eff2;

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
    //--- by default, no grid treatment
    return function_->f( (double*)&x[0], function_->dim, (void*)input_params_ );
  }

  double
  Integrator::uniform() const
  {
    return gsl_rng_uniform( rng_.get() );
  }
}

