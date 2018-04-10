#include "CepGen/Core/Integrator.h"
#include "CepGen/Core/ThreadWorker.h"
#include "CepGen/Core/utils.h"
#include "CepGen/Core/Exception.h"

#include "CepGen/Parameters.h"
#include "CepGen/Hadronisers/GenericHadroniser.h"

#include <thread>
#include <math.h>

#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_monte_miser.h>

namespace CepGen
{
  Integrator::Integrator( unsigned int ndim, double integrand( double*, size_t, void* ), Parameters* params ) :
    input_params_( params ),
    function_( new gsl_monte_function{ integrand, ndim, (void*)input_params_ } )
  {
    //--- initialise the random number generator
    gsl_rng_env_setup();
    rng_ = std::shared_ptr<gsl_rng>( gsl_rng_alloc( gsl_rng_default ), gsl_rng_free );
    unsigned long seed = ( input_params_->integrator.seed > 0 )
      ? input_params_->integrator.seed
      : time( nullptr ); // seed with time
    gsl_rng_set( rng_.get(), seed );

    input_params_->integrator.vegas.ostream = stderr; // redirect all debugging information to the error stream
    input_params_->integrator.vegas.iterations = 10;

    Debugging( Form( "Number of integration dimensions: %d,\n\t"
                     "Number of iterations [VEGAS]:     %d,\n\t"
                     "Number of function calls:         %d.",
                     function_->dim,
                     input_params_->integrator.vegas.iterations,
                     input_params_->integrator.ncvg ) );
  }

  Integrator::~Integrator()
  {}

  int
  Integrator::integrate( double& result, double& abserr )
  {
    int res = -1;
    gsl_monte_plain_state* pln_state;
    gsl_monte_vegas_state* veg_state;
    gsl_monte_miser_state* mis_state;
    const Integrator::Type algorithm = input_params_->integrator.type;

    //--- integration bounds
    std::vector<double> x_low( function_->dim, 0. ), x_up( function_->dim, 1. );

    //--- prepare integrator
    if      ( algorithm == Plain ) pln_state = gsl_monte_plain_alloc( function_->dim );
    else if ( algorithm == Vegas ) veg_state = gsl_monte_vegas_alloc( function_->dim );
    else if ( algorithm == MISER ) mis_state = gsl_monte_miser_alloc( function_->dim );

    if ( algorithm == Plain )
      res = gsl_monte_plain_integrate( function_.get(),
        &x_low[0], &x_up[0],
        function_->dim, input_params_->integrator.ncvg,
        rng_.get(), pln_state,
        &result, &abserr );
    else if ( algorithm == Vegas ) {
      gsl_monte_vegas_params_set( veg_state, &input_params_->integrator.vegas );
      //----- Vegas warmup (prepare the grid)
      if ( !grid.grid_prepared ) {
        res = gsl_monte_vegas_integrate( function_.get(),
          &x_low[0], &x_up[0],
          function_->dim, 25000,
          rng_.get(), veg_state,
          &result, &abserr );
        grid.grid_prepared = true;
      }
      Information( "Finished the Vegas warm-up." );
      //----- integration
      unsigned short it_chisq = 0;
      do {
        res = gsl_monte_vegas_integrate( function_.get(),
          &x_low[0], &x_up[0],
          function_->dim, 0.2 * input_params_->integrator.ncvg,
          rng_.get(), veg_state,
          &result, &abserr );
        PrintMessage( Form( "\t>> at call %2d: average = %10.6f   "
                            "sigma = %10.6f   chi2 = %4.3f.",
                            it_chisq+1, result, abserr,
                            gsl_monte_vegas_chisq( veg_state ) ) );
        it_chisq++;
      } while ( fabs( gsl_monte_vegas_chisq( veg_state )-1. ) > 0.5 );
    }
    //----- integration
    else if ( algorithm == MISER ) {
      gsl_monte_miser_params_set( mis_state, &input_params_->integrator.miser );
      res = gsl_monte_miser_integrate( function_.get(),
        &x_low[0], &x_up[0],
        function_->dim, input_params_->integrator.ncvg,
        rng_.get(), mis_state,
        &result, &abserr );
    }

    //--- clean integrator
    if      ( algorithm == Plain ) gsl_monte_plain_free( pln_state );
    else if ( algorithm == Vegas ) gsl_monte_vegas_free( veg_state );
    else if ( algorithm == MISER ) gsl_monte_miser_free( mis_state );

    if ( input_params_->hadroniser() )
      input_params_->hadroniser()->setCrossSection( result, abserr );

    return res;
  }

  unsigned short
  Integrator::dimensions() const
  {
    if ( !function_ )
      return 0;
    return function_->dim;
  }

  void
  Integrator::generate( unsigned long num_events, std::function<void( const Event&, unsigned long )> callback )
  {
    if ( !grid.gen_prepared )
      setGen();

    if ( input_params_->generation.num_threads > 1 )
      Information( Form( "Will generate events using %d threads", input_params_->generation.num_threads ) );

    // define the threads and workers
    std::vector<std::thread> threads;
    std::vector<std::shared_ptr<ThreadWorker> > workers;
    for ( unsigned int i = 0; i < input_params_->generation.num_threads; ++i ) {
      std::shared_ptr<ThreadWorker> worker( new ThreadWorker( &mutex_, rng_.get(), function_.get(), &grid, callback ) );
      workers.emplace_back( worker );
      threads.emplace_back( &ThreadWorker::generate, worker.get() );
    }
    // launch the multi-threaded events generation
    for ( auto& thread : threads )
      if ( thread.joinable() )
        thread.join();
  }

  void
  Integrator::setGen()
  {
    Information( Form( "Preparing the grid (%d points) for the generation of unweighted events.", input_params_->integrator.npoints ) );

    grid.max = pow( grid.mbin_, function_->dim );
    const double inv_npoin = 1./input_params_->integrator.npoints;

    if ( function_->dim > grid.max_dimensions_ )
      FatalError( Form( "Number of dimensions to integrate exceeds the maximum number, %d", grid.max_dimensions_ ) );

    grid.f_max = std::vector<double>( grid.max, 0. );
    grid.n = std::vector<int>( function_->dim, 0 );

    std::vector<double> x( function_->dim, 0. );

    input_params_->generation.ngen = 0;
    input_params_->setStorage( false );

    // ...
    double sum = 0., sum2 = 0., sum2p = 0.;

    //--- main loop
    for ( unsigned int i = 0; i < grid.max; ++i ) {
      int jj = i;
      for ( unsigned int j = 0; j < function_->dim; ++j ) {
        int jjj = jj*grid.inv_mbin_;
        grid.n[j] = jj-jjj*grid.mbin_;
        jj = jjj;
      }
      double fsum = 0., fsum2 = 0.;
      for ( unsigned int j = 0; j < input_params_->integrator.npoints; ++j ) {
        for ( unsigned int k = 0; k < function_->dim; ++k )
          x[k] = ( gsl_rng_uniform( rng_.get() ) + grid.n[k] ) * grid.inv_mbin_;
        const double z = function_->f( (double*)&x[0], function_->dim, (void*)input_params_ );
        grid.f_max[i] = std::max( grid.f_max[i], z );
        fsum += z;
        fsum2 += z*z;
      }
      const double av = fsum*inv_npoin, av2 = fsum2*inv_npoin, sig2 = av2 - av*av;
      sum += av;
      sum2 += av2;
      sum2p += sig2;
      grid.f_max_global = std::max( grid.f_max_global, grid.f_max[i] );

      if ( Logger::get().level >= Logger::DebugInsideLoop ) {
        const double sig = sqrt( sig2 );
        const double eff = ( grid.f_max[i] != 0. ) ? grid.f_max[i]/av : 1.e4;
        std::ostringstream os;
        for ( unsigned int j = 0; j < function_->dim; ++j )
          os << grid.n[j] << ( j != function_->dim-1 ? ", " : "" );
        DebuggingInsideLoop( Form( "In iteration #%d:\n\t"
                                   "av   = %f\n\t"
                                   "sig  = %f\n\t"
                                   "fmax = %f\n\t"
                                   "eff  = %f\n\t"
                                   "n = (%s)",
                                   i, av, sig, grid.f_max[i], eff, os.str().c_str() ) );
      }
    } // end of main loop

    const double inv_max = 1./grid.max;
    sum *= inv_max;
    sum2 *= inv_max;
    sum2p *= inv_max;

    if ( Logger::get().level >= Logger::Debug ) {
      const double sig = sqrt( sum2-sum*sum ), sigp = sqrt( sum2p );

      double eff1 = 0.;
      for ( unsigned int i = 0; i < grid.max; ++i )
        eff1 += ( grid.f_max[i] / ( grid.max*sum ) );
      const double eff2 = grid.f_max_global/sum;

      Debugging( Form( "Average function value     = sum   = %g\n\t"
                       "Average function value**2  = sum2  = %g\n\t"
                       "Overall standard deviation = sig   = %g\n\t"
                       "Average standard deviation = sigp  = %g\n\t"
                       "Maximum function value     = f_max = %g\n\t"
                       "Average inefficiency       = eff1  = %g\n\t"
                       "Overall inefficiency       = eff2  = %g\n\t",
                       sum, sum2, sig, sigp, grid.f_max_global, eff1, eff2 ) );
    }
    grid.gen_prepared = true;
    Information( "Grid prepared! Now launching the production." );
  }

  std::ostream&
  operator<<( std::ostream& os, const Integrator::Type& type )
  {
    switch ( type ) {
      case Integrator::Plain: return os << "Plain";
      case Integrator::Vegas: return os << "Vegas";
      case Integrator::MISER: return os << "MISER";
    }
    return os;
  }

  //------------------------------------------------------------------------------------------------

  GridParameters::GridParameters() :
    grid_prepared( false ), gen_prepared( false ),
    f_max_global( 0. )
  {}
}

