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
    ps_bin_( 0 ), function_( function ), grid_( grid ),
    grid_correc_( 0. ), grid_correc2_( 0. ),
    grid_f_max2_( 0. ), grid_f_max_diff_( 0. ), grid_f_max_old_( 0. ),
    mutex_( mutex ), callback_( callback )
  {
    if ( !function )
      throw FatalError( "ThreadWorker" ) << "Invalid integration function passed!";

    rng_ = std::shared_ptr<gsl_rng>( gsl_rng_clone( rng ), gsl_rng_free );
    grid_nm_.reserve( grid_->max );
    // retrieve standard parameters
    global_params_ = static_cast<Parameters*>( function->params );
    // clone the process for this thread
    process_ = global_params_->processClone();
    // copy the standard parameters and feed the cloned process
    local_params_ = new Parameters( *static_cast<const Parameters*>( global_params_ ) );
    local_params_->setProcess( process_.get() );
  }

  bool
  ThreadWorker::generate()
  {
    if ( !grid_->gen_prepared )
      throw FatalError( "ThreadWorker" ) << "Generation not prepared!";

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

    std::vector<double> x( function_->dim, 0. );

    //--- correction cycles

    if ( ps_bin_ != 0 ) {
      bool has_correction = false;
      while ( !correctionCycle( x, has_correction ) ) {/*std::cout <<"correctioncycle"<<std::endl;*/}
      if ( has_correction )
        return storeEvent( x );
    }

    double weight = 0., y = -1.;

    //--- normal generation cycle

    //----- select a Integrator bin and reject if fmax is too small
    do {
      do {
        // ...
        ps_bin_ = uniform() * grid_->max;
        y = uniform() * grid_->f_max_global;
        grid_nm_[ps_bin_] += 1;
      } while ( y > grid_->f_max[ps_bin_] );
      // shoot a point x in this bin
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
//std::cout << weight << "|" << y << std::endl;
    } while ( y > weight );

//std::cout << __PRETTY_FUNCTION__<<"|"<< weight<< std::endl;

    if ( weight <= grid_->f_max[ps_bin_] )
      ps_bin_ = 0;
    // Init correction cycle if weight is higher than fmax or ffmax
    else if ( weight <= grid_->f_max_global ) {
      grid_f_max_old_ = grid_->f_max[ps_bin_];
      grid_->f_max[ps_bin_] = weight;
      grid_f_max_diff_ = weight-grid_f_max_old_;
      grid_correc_ = ( grid_nm_[ps_bin_] - 1. ) * grid_f_max_diff_ / grid_->f_max_global - 1.;
    }
    else {
      grid_f_max_old_ = grid_->f_max[ps_bin_];
      grid_->f_max[ps_bin_] = weight;
      grid_f_max_diff_ = weight-grid_f_max_old_;
      grid_->f_max_global = weight;
      grid_correc_ = ( grid_nm_[ps_bin_] - 1. ) * grid_f_max_diff_ / grid_->f_max_global * weight / grid_->f_max_global - 1.;
    }

    DebuggingInsideLoop( "ThreadWorker" )
      << "Correction applied: " << grid_correc_ << ", phase space bin = " << ps_bin_;

    // Return with an accepted event
    if ( weight > 0. )
      return storeEvent( x );
    return false;
  }

  bool
  ThreadWorker::correctionCycle( std::vector<double>& x, bool& has_correction )
  {
//    std::cout << __PRETTY_FUNCTION__ << std::endl;
    DebuggingInsideLoop( "ThreadWorker" )
      << "Correction cycles are started.\n\t"
      << "bin = " << ps_bin_ << "\t"
      << "correc = " << grid_correc_ << "\t"
      << "corre2 = " << grid_correc2_ << ".";

    if ( grid_correc_ >= 1. )
      grid_correc_ -= 1.;
    if ( uniform() < grid_correc_ ) {
      grid_correc_ = -1.;
      std::vector<double> xtmp( function_->dim );
      // Select x values in phase space bin
      for ( unsigned int k = 0; k < function_->dim; ++k )
        xtmp[k] = ( uniform() + grid_->n[k] ) * grid_->inv_mbin_;
      // Compute weight for x value
      const double weight = eval( xtmp );
      // Parameter for correction of correction
      if ( weight > grid_->f_max[ps_bin_] ) {
        if ( weight > grid_f_max2_ )
          grid_f_max2_ = weight;
        grid_correc2_ -= 1.;
        grid_correc_ += 1.;
      }
      // Accept event
      if ( weight >= grid_f_max_diff_*uniform() + grid_f_max_old_ ) { // FIXME!!!!
//        InError("Accepting event!!!");
        //return storeEvent(x);
        x = xtmp;
        has_correction = true;
        return true;
      }
      return false;
    }
    // Correction if too big weight is found while correction
    // (All your bases are belong to us...)
    if ( grid_f_max2_ > grid_->f_max[ps_bin_] ) {
      grid_f_max_old_ = grid_->f_max[ps_bin_];
      grid_->f_max[ps_bin_] = grid_f_max2_;
      grid_f_max_diff_ = grid_f_max2_-grid_f_max_old_;
      const double correc_tmp = ( grid_nm_[ps_bin_]-1. ) * grid_f_max_diff_ / grid_->f_max_global;
      if ( grid_f_max2_ < grid_->f_max_global )
        grid_correc_ = correc_tmp - grid_correc2_;
      else {
        grid_->f_max_global = grid_f_max2_;
        grid_correc_ = correc_tmp * grid_f_max2_ / grid_->f_max_global - grid_correc2_;
      }
      grid_correc2_ = 0.;
      grid_f_max2_ = 0.;
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
    return gsl_rng_uniform( rng_.get() );
  }

  bool
  ThreadWorker::storeEvent( const std::vector<double>& x )
  {
//    std::cout << "------------store!!!" << std::endl;
    std::lock_guard<std::mutex> guard( *mutex_ );
    //mutex_->lock();
    local_params_->setStorage( true );
    const double weight = eval( x );
    local_params_->setStorage( false );
    //mutex_->unlock();

    if ( weight <= 0. )
      return false;

    if ( global_params_->generation.ngen % global_params_->generation.gen_print_every == 0 ) {
      Information( "ThreadWorker" )
        << "[thread 0x" << std::hex << std::hash<std::thread::id>()( std::this_thread::get_id() ) << std::dec
        << "] Generated events: " << global_params_->generation.ngen;
      local_params_->process()->last_event->dump();
    }
    if ( callback_ )
      callback_( *local_params_->process()->last_event, global_params_->generation.ngen );

    global_params_->generation.ngen += 1;
    return true;
  }
}

