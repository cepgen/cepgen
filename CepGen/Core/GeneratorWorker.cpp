#include "CepGen/Core/GeneratorWorker.h"

#include "CepGen/Integration/Integrator.h"
#include "CepGen/Integration/Integrand.h"
#include "CepGen/Integration/GridParameters.h"

#include "CepGen/Core/EventModifier.h"
#include "CepGen/Core/ExportModule.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/Exception.h"

#include "CepGen/Parameters.h"
#include "CepGen/Processes/Process.h"
#include "CepGen/Event/Event.h"

#include "CepGen/Utils/String.h"
#include "CepGen/Utils/ProgressBar.h"

namespace cepgen
{
  GeneratorWorker::GeneratorWorker( Parameters* params ) :
    integrand_( new Integrand( params ) ), integrator_( nullptr ), params_( params ),
    ps_bin_( INVALID_BIN )
  {
    CG_INFO( "GeneratorWorker" )
      << "New generator worker initialised for integration/event generation.\n\t"
      << "Parameters at " << (void*)params_ << ".";
  }

  void
  GeneratorWorker::setIntegrator( const Integrator* integr )
  {
    integrator_ = integr;
    grid_.reset( new GridParameters( integrator_->size() ) );
    CG_DEBUG( "GeneratorWorker:integrator" )
      << "Dim-" << integrator_->size() << " " << integrator_->name() << " integrator "
      << "set for dim-" << grid_->n( 0 ).size() << " grid.";
  }

  //-----------------------------------------------------------------------------------------------
  // events generation part
  //-----------------------------------------------------------------------------------------------

  void
  GeneratorWorker::generate( size_t num_events, Event::callback callback )
  {
    if ( !params_ )
      throw CG_FATAL( "GeneratorWorker:generate" ) << "No steering parameters specified!";

    if ( num_events < 1 )
      num_events = params_->generation().maxgen;

    try {
      while ( params_->numGeneratedEvents() < num_events )
        next( callback );
    } catch ( const Exception& e ) {
      CG_WARNING( "GeneratorWorker:generate" ) << "Generation ended with exception.";
      e.dump();
      return;
    }
  }

  bool
  GeneratorWorker::next( Event::callback callback )
  {
    if ( !integrator_ )
      throw CG_FATAL( "GeneratorWorker:generate" ) << "No integrator object handled!";

    //--- a few checks on the grid
    if ( !grid_ )
      throw CG_FATAL( "GeneratorWorker:generate" ) << "No grid object handled!";
    if ( !grid_->gen_prepared )
      computeGenerationParameters();

    std::vector<double> xtmp( integrator_->size() );

    //--- correction cycles

    if ( ps_bin_ != INVALID_BIN ) {
      bool has_correction = false;
      while ( !correctionCycle( xtmp, has_correction ) ) {}
      if ( has_correction )
        return storeEvent( xtmp, callback );
    }

    double weight = 0.;

    //--- normal generation cycle

    while ( true ) {
      double y = -1.;
      //----- select a and reject if fmax is too small
      while ( true ) {
        // ...
        ps_bin_ = integrator_->uniform() * grid_->size();
        y = integrator_->uniform() * grid_->globalMax();
        grid_->setTrial( ps_bin_ );
        if ( y <= grid_->maxValue( ps_bin_ ) )
          break;
      }
      // shoot a point x in this bin
      grid_->shoot( integrator_, ps_bin_, xtmp );
      // get weight for selected x value
      weight = integrator_->eval( xtmp );
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

      CG_DEBUG("GeneratorWorker:generateOne")
        << "Correction " << grid_->correc << " will be applied "
        << "for phase space bin " << ps_bin_ << ".";
    }

    // return with an accepted event
    return storeEvent( xtmp, callback );
  }

  bool
  GeneratorWorker::correctionCycle( std::vector<double>& x, bool& has_correction )
  {
    CG_DEBUG_LOOP( "GeneratorWorker:correction" )
      << "Correction cycles are started.\n\t"
      << "bin = " << ps_bin_ << "\t"
      << "correc = " << grid_->correc << "\t"
      << "corre2 = " << grid_->correc2 << ".";

    if ( grid_->correc >= 1. )
      grid_->correc -= 1.;

    if ( integrator_->uniform() < grid_->correc ) {
      grid_->correc = -1.;
      std::vector<double> xtmp( integrator_->size() );
      // Select x values in phase space bin
      grid_->shoot( integrator_, ps_bin_, xtmp );
      const double weight = integrator_->eval( xtmp );
      // Parameter for correction of correction
      if ( weight > grid_->maxValue( ps_bin_ ) ) {
        grid_->f_max2 = std::max( grid_->f_max2, weight );
        grid_->correc += 1.;
        grid_->correc2 -= 1.;
      }
      // Accept event
      if ( weight >= grid_->f_max_diff*integrator_->uniform() + grid_->f_max_old ) {
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
  GeneratorWorker::storeEvent( const std::vector<double>& x, Event::callback callback )
  {
    //--- start by computing the matrix element for that point
    const double weight = integrator_->eval( x );

    //--- reject if unphysical
    if ( weight <= 0. )
      return false;

    const auto ngen = params_->numGeneratedEvents();
    if ( integrand_->process().hasEvent() ) {
      auto& event = integrand_->process().event();
      if ( ( ngen+1 ) % params_->generation().gen_print_every == 0 )
        CG_INFO( "GeneratorWorker:store" )
          << "Generated events: " << ngen+1 << "\n" << event;
      if ( callback )
        callback( event, ngen );
      for ( auto& mod : params_->outputModulesSequence() )
        *mod << event;
      const_cast<Parameters*>( params_ )->addGenerationTime( event.time_total );
    }
    return true;
  }

  //-----------------------------------------------------------------------------------------------
  // initial preparation run before the generation of unweighted events
  //-----------------------------------------------------------------------------------------------

  void
  GeneratorWorker::computeGenerationParameters()
  {
    if ( !params_ )
      throw CG_FATAL( "GeneratorWorker:setGen" ) << "No steering parameters specified!";

    integrand_->setStorage( false );

    CG_INFO( "GeneratorWorker:setGen" )
      << "Preparing the grid (" << params_->generation().num_points << " points/bin) "
      << "for the generation of unweighted events.";

    const double inv_num_points = 1./params_->generation().num_points;
    std::vector<double> point_coord( integrator_->size(), 0. );
    if ( point_coord.size() < grid_->n( 0 ).size() )
      throw CG_FATAL( "GridParameters:shoot" )
        << "Coordinates vector multiplicity is insufficient!";

    // ...
    double sum = 0., sum2 = 0., sum2p = 0.;

    utils::ProgressBar prog_bar( grid_->size(), 5 );

    //--- main loop
    for ( unsigned int i = 0; i < grid_->size(); ++i ) {
      double fsum = 0., fsum2 = 0.;
      for ( unsigned int j = 0; j < params_->generation().num_points; ++j ) {
        grid_->shoot( integrator_, i, point_coord );
        const double weight = integrator_->eval( point_coord );
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
      if ( CG_LOG_MATCH( "GeneratorWorker:setGen", debugInsideLoop ) ) {
        const double sig = sqrt( sig2 );
        const double eff = ( grid_->maxValue( i ) != 0. )
          ? av/grid_->maxValue( i )
          : 0.;
        CG_DEBUG_LOOP( "GeneratorWorker:setGen" )
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

    CG_DEBUG( "GeneratorWorker:setGen" )
      << "Average function value         = " << sum << "\n\t"
      << "Average squared function value = " << sum2 << "\n\t"
      << "Overall standard deviation     = " << sig << "\n\t"
      << "Average standard deviation     = " << sigp << "\n\t"
      << "Maximum function value         = " << grid_->globalMax() << "\n\t"
      << "Average inefficiency           = " << eff1 << "\n\t"
      << "Overall inefficiency           = " << eff2;

    grid_->gen_prepared = true;
    //--- from now on events will be stored
    integrand_->setStorage( true );

    CG_INFO( "GeneratorWorker:setGen" )
      << "Grid prepared! Now launching the production.";
  }
}

