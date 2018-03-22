#include "CepGen/Core/ThreadWorker.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"

#include "CepGen/Event/Event.h"
#include "CepGen/Processes/GenericProcess.h"

#include "CepGen/Parameters.h"

#include <thread>
#include <math.h>

namespace CepGen
{
  extern int gSignal;

  ThreadWorker::ThreadWorker( std::mutex* mutex,
                              gsl_rng* rng, gsl_monte_function* function,
                              GridParameters* grid,
                              std::function<void( const Event&, unsigned long )>& callback ) :
    ps_bin_( 0 ), rng_( rng ), function_( function ), grid_( grid ),
    mutex_( mutex ), callback_( callback )
  {
    if ( !function )
      throw Exception( __PRETTY_FUNCTION__, "Invalid integration function passed!", FatalError );

    global_params_ = static_cast<Parameters*>( function->params );
    process_ = global_params_->processClone();
    local_params_ = new Parameters( *static_cast<const Parameters*>( global_params_ ) );
    local_params_->setProcess( process_.get() );
  }

  bool
  ThreadWorker::generate()
  {
    if ( !grid_->gen_prepared )
      throw Exception( __PRETTY_FUNCTION__, "Generation not prepared!", FatalError );

    while ( true ) {
      if ( !next() )
        continue;
      if ( gSignal != 0 )
        return false;
      if ( global_params_->generation.ngen >= global_params_->generation.maxgen )
        return true;
    }
    return false;
  }

  bool
  ThreadWorker::next()
  {
    if ( global_params_->generation.ngen >= global_params_->generation.maxgen )
      return true;

    const unsigned int max = pow( grid_->mbin_, function_->dim );

    std::vector<double> x( function_->dim, 0. );

    //--- correction cycles
    
    if ( ps_bin_ != 0 ) {
      bool has_correction = false;
      while ( !correctionCycle( x, has_correction ) ) {}
      if ( has_correction )
        return storeEvent( x );
    }

    double weight = 0., y = -1.;

    //--- normal generation cycle

    //----- select a Integrator bin and reject if fmax is too small
    do {
      do {
        // ...
        ps_bin_ = uniform() * max;
        y = uniform() * grid_->f_max_global;
        grid_->nm[ps_bin_] += 1;
      } while ( y > grid_->f_max[ps_bin_] );
      // Select x values in this Integrator bin
      int jj = ps_bin_;
      for ( unsigned int i = 0; i < function_->dim; ++i ) {
        int jjj = jj * grid_->inv_mbin_;
        grid_->n[i] = jj - jjj * grid_->mbin_;
        x[i] = ( uniform() + grid_->n[i] ) * grid_->inv_mbin_;
        jj = jjj;
      }

      // Get weight for selected x value
      weight = eval( x );
      if ( weight <= 0. )
        continue;
    } while ( y > weight );

    if ( weight <= grid_->f_max[ps_bin_] )
      ps_bin_ = 0;
    // Init correction cycle if weight is higher than fmax or ffmax
    else if ( weight <= grid_->f_max_global ) {
      grid_->f_max_old = grid_->f_max[ps_bin_];
      grid_->f_max[ps_bin_] = weight;
      grid_->f_max_diff = weight-grid_->f_max_old;
      grid_->correc = ( grid_->nm[ps_bin_] - 1. ) * grid_->f_max_diff / grid_->f_max_global - 1.;
    }
    else {
      grid_->f_max_old = grid_->f_max[ps_bin_];
      grid_->f_max[ps_bin_] = weight;
      grid_->f_max_diff = weight-grid_->f_max_old;
      grid_->f_max_global = weight;
      grid_->correc = ( grid_->nm[ps_bin_] - 1. ) * grid_->f_max_diff / grid_->f_max_global * weight / grid_->f_max_global - 1.;
    }

    DebuggingInsideLoop( Form( "Correction applied: %f, phase space bin = %d", grid_->correc, ps_bin_ ) );

    // Return with an accepted event
    if ( weight > 0. )
      return storeEvent( x );
    return false;
  }

  bool
  ThreadWorker::correctionCycle( std::vector<double>& x, bool& has_correction )
  {
    DebuggingInsideLoop( Form( "Correction cycles are started.\n\t"
                               "bin = %d\t"
                               "correc = %g\t"
                               "corre2 = %g.", ps_bin_, grid_->correc, grid_->correc2 ) );

    if ( grid_->correc >= 1. )
      grid_->correc -= 1.;
    if ( uniform() < grid_->correc ) {
      grid_->correc = -1.;
      std::vector<double> xtmp;
      // Select x values in phase space bin
      for ( unsigned int k = 0; k < function_->dim; ++k )
        xtmp.emplace_back( ( uniform() + grid_->n[k] ) * grid_->inv_mbin_ );
      // Compute weight for x value
      const double weight = eval( xtmp );
      // Parameter for correction of correction
      if ( weight > grid_->f_max[ps_bin_] ) {
        if ( weight > grid_->f_max2 )
          grid_->f_max2 = weight;
        grid_->correc2 -= 1.;
        grid_->correc += 1.;
      }
      // Accept event
      if ( weight >= grid_->f_max_diff*uniform() + grid_->f_max_old ) { // FIXME!!!!
        //Error("Accepting event!!!");
        //return storeEvent(x);
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
      const double correc_tmp = ( grid_->nm[ps_bin_] - 1. ) * grid_->f_max_diff / grid_->f_max_global;
      if ( grid_->f_max2 < grid_->f_max_global )
        grid_->correc = correc_tmp - grid_->correc2;
      else {
        grid_->f_max_global = grid_->f_max2;
        grid_->correc = correc_tmp * grid_->f_max2 / grid_->f_max_global - grid_->correc2;
      }
      grid_->correc2 = 0.;
      grid_->f_max2 = 0.;
      return false;
    }
    return true;
  }

  double
  ThreadWorker::eval( const std::vector<double>& x )
  {
    return function_->f( (double*)&x[0], function_->dim, (void*)local_params_ );
  }

  double
  ThreadWorker::uniform() const
  {
    return gsl_rng_uniform( rng_ );
  }

  bool
  ThreadWorker::storeEvent( const std::vector<double>& x )
  {
    std::lock_guard<std::mutex> guard( *mutex_ );
    local_params_->setStorage( true );
    const double weight = eval( x );
    local_params_->setStorage( false );

    if ( weight <= 0. )
      return false;

    global_params_->generation.ngen += 1;

    if ( global_params_->generation.ngen % global_params_->generation.gen_print_every == 0 ) {
      //mutex_->lock();
      Information( Form( "[thread 0x%zx] Generated events: %d",
                         std::hash<std::thread::id>()( std::this_thread::get_id() ),
                         global_params_->generation.ngen ) );
      local_params_->process()->last_event->dump();
      //mutex_->unlock();
    }

    if ( callback_ ) {
      callback_( *local_params_->process()->last_event, global_params_->generation.ngen );
    }

    return true;
  }
}
