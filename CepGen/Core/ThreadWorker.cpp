#include "CepGen/Core/ThreadWorker.h"
#include "CepGen/Core/GridParameters.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"

#include "CepGen/Event/Event.h"
#include "CepGen/Processes/GenericProcess.h"
#include "CepGen/Processes/GenericKTProcess.h" //FIXME

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
    ps_bin_( 0 ), rng_( rng ), function_( function ),
    grid_( grid ), grid_num_( grid_->max, 0 ),
    grid_correc_( 0. ), grid_correc2_( 0. ),
    grid_f_max2_( 0. ), /*grid_f_max_diff_( 0. ), */grid_f_max_old_( 0. ),
    global_params_( nullptr ),
    mutex_( mutex ), callback_( callback )
  {
    if ( !function )
      throw CG_FATAL( "ThreadWorker" ) << "Invalid integration function passed!";

    // retrieve standard parameters
    global_params_ = static_cast<Parameters*>( function->params );
    // copy the standard parameters and feed the cloned process
    local_params_ = std::unique_ptr<Parameters>( new Parameters( *const_cast<const Parameters*>( global_params_ ) ) );
    // clone the process for this thread
    local_params_->cloneProcess( global_params_->process() );
  }

  bool
  ThreadWorker::generate( unsigned long max_gen )
  {
    if ( !grid_->gen_prepared )
      throw CG_FATAL( "ThreadWorker" ) << "Generation not prepared!";

    grid_correc_ = 0.;
    while ( true ) {
      // only keep physical events
      if ( !next() )
        continue;
      // check if the user interrupted the generation
      if ( gSignal != 0 )
        return false;
      // check if we generated enough events for this thread
      if ( max_gen > 0 && local_params_->generation.ngen >= max_gen )
        return true;
      // check if we generated enough events for the full run [all threads]
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
      while ( !correctionCycle( x, has_correction ) ) {}
      if ( has_correction )
        return storeEvent( x );
    }

    //--- normal generation cycle

    double weight = 0.;

    while ( gSignal == 0 ) {
      double y = -1.;
      //--- select a bin and reject if fmax is too small
      while ( gSignal == 0 ) {
        ps_bin_ = uniform() * grid_->max;
        grid_num_[ps_bin_]++;
        y = uniform() * grid_->f_max_global;
        if ( y <= grid_->f_max[ps_bin_] )
          break;
      }
      //std::cout << "---> " << y << "::" << ps_bin_ << std::endl;
      // shoot a point x in this bin
      std::vector<unsigned short> grid_n = grid_->n_map.at( ps_bin_ );
      for ( unsigned int i = 0; i < function_->dim; ++i )
        x[i] = ( uniform()+grid_n[i] ) * GridParameters::inv_mbin_;
      // get weight for selected x value
      weight = eval( x );
//std::cout << weight << std::endl;
      if ( weight >= y )
        break;
//std::cout << weight << "|" << y << std::endl;
    }
    if ( gSignal != 0 )
      return false;

//std::cout << __PRETTY_FUNCTION__<<"|"<< weight<< std::endl;

    if ( weight < grid_->f_max[ps_bin_] )
      ps_bin_ = 0;
    // init correction cycle if weight is higher than local or global maximum
    else if ( weight <= grid_->f_max_global ) {
      grid_f_max_old_ = grid_->f_max[ps_bin_];
      grid_->f_max[ps_bin_] = weight;
      grid_->f_max_diff = weight-grid_f_max_old_;
      grid_correc_ = ( grid_num_[ps_bin_]-1. ) * grid_->f_max_diff / grid_->f_max_global - 1.;
    }
    else {
      // a new function global maximum has been found!
      // grid correction is hence needed
      grid_f_max_old_ = grid_->f_max[ps_bin_];
      grid_->f_max[ps_bin_] = weight;
      grid_->f_max_diff = weight-grid_f_max_old_;
      grid_->f_max_global = weight;
      grid_correc_ = ( grid_num_[ps_bin_]-1. ) * grid_->f_max_diff / grid_->f_max_global * weight / grid_->f_max_global - 1.;
    }

    CG_DEBUG_LOOP( "ThreadWorker:next" )
      << "Correction " << grid_correc_ << " will be applied for phase space bin " << ps_bin_ << ".";

    // return with an accepted event
    if ( weight > 0. )
      return storeEvent( x );
    return false;
  }

  bool
  ThreadWorker::correctionCycle( std::vector<double>& x, bool& has_correction )
  {
    CG_DEBUG_LOOP( "ThreadWorker:correction" )
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
      const std::vector<unsigned short> grid_n = grid_->n_map.at( ps_bin_ );
      for ( unsigned int k = 0; k < function_->dim; ++k )
        xtmp[k] = ( uniform() + grid_n[k] ) * GridParameters::inv_mbin_;
      const double weight = eval( xtmp );
      // Parameter for correction of correction
      if ( weight > grid_->f_max[ps_bin_] ) {
        grid_f_max2_ = std::max( grid_f_max2_, weight );
        grid_correc_ += 1.;
        grid_correc2_ -= 1.;
      }
      // Accept event
      if ( weight >= grid_->f_max_diff*uniform() + grid_f_max_old_ ) {
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
      grid_->f_max_diff = grid_f_max2_-grid_f_max_old_;
      const double correc_tmp = ( grid_num_[ps_bin_]-1. ) * grid_->f_max_diff / grid_->f_max_global;
      if ( grid_f_max2_ < grid_->f_max_global )
        grid_correc_ = correc_tmp;
      else {
        grid_->f_max_global = grid_f_max2_;
        grid_correc_ = correc_tmp * grid_f_max2_ / grid_->f_max_global;
      }
      grid_correc_ -= grid_correc2_;
      grid_correc2_ = 0.;
      grid_f_max2_ = 0.;
      return false;
    }
    return true;
  }

  bool
  ThreadWorker::storeEvent( const std::vector<double>& x )
  {
//    std::cout << "------------store!!! " << global_params_->generation.ngen << std::endl;
    local_params_->setStorage( true );
    const double weight = eval( x );
    local_params_->setStorage( false );

    if ( weight <= 0. )
      return false;

    {
      std::lock_guard<std::mutex> guard( *mutex_ );
      if ( global_params_->generation.ngen % global_params_->generation.gen_print_every == 0 ) {
        CG_INFO( "ThreadWorker:store" )
          << "[thread 0x" << std::hex << std::hash<std::thread::id>()( std::this_thread::get_id() ) << std::dec
          << "] Generated events: " << global_params_->generation.ngen;
        local_params_->process()->last_event->dump();
      }
      global_params_->process()->last_event = local_params_->process()->last_event;

      local_params_->generation.ngen += 1;
      global_params_->generation.ngen += 1;

      if ( callback_ )
        callback_( *local_params_->process()->last_event, global_params_->generation.ngen );
    }
    return true;
  }

  //-----------------------------------------------------------------------------------------------
  // Helper methods
  //-----------------------------------------------------------------------------------------------

  double
  ThreadWorker::eval( const std::vector<double>& x )
  {
    return function_->f( (double*)&x[0], function_->dim, (void*)local_params_.get() );
  }

  double
  ThreadWorker::uniform() const
  {
    return gsl_rng_uniform( rng_ );
  }
}

