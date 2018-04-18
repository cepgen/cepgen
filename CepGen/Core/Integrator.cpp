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
    rng_ = std::shared_ptr<gsl_rng>( gsl_rng_alloc( input_params_->integrator.rng_engine ), gsl_rng_free );
    unsigned long seed = ( input_params_->integrator.rng_seed > 0 )
      ? input_params_->integrator.rng_seed
      : time( nullptr ); // seed with time
    gsl_rng_set( rng_.get(), seed );

    input_params_->integrator.vegas.ostream = stderr; // redirect all debugging information to the error stream
    input_params_->integrator.vegas.iterations = 10;

    CG_DEBUG( "Integrator:build" )
      << "Number of integration dimensions: " << function_->dim << ",\n\t"
      << "Number of iterations [VEGAS]:     " << input_params_->integrator.vegas.iterations << ",\n\t"
      << "Number of function calls:         " << input_params_->integrator.ncvg << ",\n\t"
      << "Random numbers generator:         " << gsl_rng_name( rng_.get() ) << ".";
  }

  Integrator::~Integrator()
  {}

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
    if ( algorithm == Plain )
      pln_state = gsl_monte_plain_alloc( function_->dim );
    else if ( algorithm == Vegas )
      veg_state = gsl_monte_vegas_alloc( function_->dim );
    else if ( algorithm == MISER )
      mis_state = gsl_monte_miser_alloc( function_->dim );

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
    else if ( algorithm == MISER ) {
      gsl_monte_miser_params_set( mis_state, &input_params_->integrator.miser );
      res = gsl_monte_miser_integrate( function_.get(),
        &x_low[0], &x_up[0],
        function_->dim, input_params_->integrator.ncvg,
        rng_.get(), mis_state,
        &result, &abserr );
    }

    //--- clean integrator
    if ( algorithm == Plain )
      gsl_monte_plain_free( pln_state );
    else if ( algorithm == Vegas )
      gsl_monte_vegas_free( veg_state );
    else if ( algorithm == MISER )
      gsl_monte_miser_free( mis_state );

    input_params_->integrator.result = result;
    input_params_->integrator.err_result = abserr;

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

  double
  Integrator::eval( const std::vector<double>& x )
  {
    return function_->f( (double*)&x[0], function_->dim, (void*)input_params_ );
  }

  void
  Integrator::generate( unsigned long num_events, std::function<void( const Event&, unsigned long )> callback )
  {
    if ( !grid.gen_prepared )
      computeGenerationParameters();

    if ( input_params_->generation.num_threads > 1 )
      CG_INFO( "Integrator:generate" )
        << "Will generate events using " << input_params_->generation.num_threads << " threads.";

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
  Integrator::computeGenerationParameters()
  {
    CG_INFO( "Integrator:setGen" )
      << "Preparing the grid (" << input_params_->integrator.npoints << " points) "
      << "for the generation of unweighted events.";

    grid.max = pow( grid.mbin_, function_->dim );
    const double inv_npoin = 1./input_params_->integrator.npoints;

    if ( function_->dim > grid.max_dimensions_ )
      throw CG_FATAL( "Integrator:setGen" )
       << "Number of dimensions to integrate exceeds the maximum number, " << grid.max_dimensions_ << ".";

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
        const double z = eval( x );
        grid.f_max[i] = std::max( grid.f_max[i], z );
        fsum += z;
        fsum2 += z*z;
      }
      const double av = fsum*inv_npoin, av2 = fsum2*inv_npoin, sig2 = av2-av*av;
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
        CG_DEBUG_LOOP( "Integrator" )
          << "In iteration #" << i << ":\n\t"
          << "av   = " << av << "\n\t"
          << "sig  = " << sig << "\n\t"
          << "fmax = " << grid.f_max[i] << "\n\t"
          << "eff  = " << eff << "\n\t"
          << "n = (" << os.str() << ")";
      }
    } // end of main loop

    const double inv_max = 1./grid.max;
    sum *= inv_max;
    sum2 *= inv_max;
    sum2p *= inv_max;

    const double sig = sqrt( sum2-sum*sum ), sigp = sqrt( sum2p );

    double eff1 = 0.;
    for ( unsigned int i = 0; i < grid.max; ++i )
      eff1 += sum*grid.max/grid.f_max[i];
    const double eff2 = sum/grid.f_max_global;

    CG_DEBUG( "Integrator:setGen" )
      << "Average function value     = sum   = " << sum << "\n\t"
      << "Average function value**2  = sum2  = " << sum2 << "\n\t"
      << "Overall standard deviation = sig   = " << sig << "\n\t"
      << "Average standard deviation = sigp  = " << sigp << "\n\t"
      << "Maximum function value     = f_max = " << grid.f_max_global << "\n\t"
      << "Average inefficiency       = eff1  = " << eff1 << "\n\t"
      << "Overall inefficiency       = eff2  = " << eff2;

    grid.gen_prepared = true;
    CG_INFO( "Integrator:setGen" ) << "Grid prepared! Now launching the production.";
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

