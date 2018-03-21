#include "Integrator.h"

#include "CepGen/Parameters.h"
#include "CepGen/Hadronisers/GenericHadroniser.h"

#include "CepGen/Core/utils.h"
#include "CepGen/Core/Exception.h"

#include "CepGen/Event/Event.h"

#include <fstream>
#include <thread>
#include <math.h>

#ifdef TIMING_ANALYSIS
extern std::vector<double> steps;
extern unsigned long num_points;
#endif

namespace CepGen
{
  Integrator::Integrator( const unsigned int dim, double f_( double*, size_t, void* ), Parameters* param ) :
    input_params_( param ),
    function_( std::unique_ptr<gsl_monte_function>( new gsl_monte_function ) )
  {
    //--- function to be integrated
    function_->f = f_;
    function_->dim = dim;
    function_->params = (void*)param;

    //--- initialise the random number generator
    gsl_rng_env_setup();
    rng_ = gsl_rng_alloc( gsl_rng_default );
    unsigned long seed = ( param->integrator.seed > 0 )
      ? param->integrator.seed
      : time( nullptr ); // seed with time
    gsl_rng_set( rng_, seed );

    input_params_->integrator.vegas.ostream = stderr; // redirect all debugging information to the error stream
    input_params_->integrator.vegas.iterations = 10;

    Debugging( Form( "Number of integration dimensions: %d,\n\t"
                     "Number of iterations [VEGAS]:     %d,\n\t"
                     "Number of function calls:         %d.",
                     dim,
                     input_params_->integrator.vegas.iterations,
                     input_params_->integrator.ncvg ) );
  }

  Integrator::~Integrator()
  {
    if ( rng_ )
      gsl_rng_free( rng_ );
  }

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
        rng_, pln_state,
        &result, &abserr );
    else if ( algorithm == Vegas ) {
      gsl_monte_vegas_params_set( veg_state, &input_params_->integrator.vegas );
      //----- Vegas warmup (prepare the grid)
      if ( !grid_.grid_prepared ) {
        res = gsl_monte_vegas_integrate( function_.get(),
          &x_low[0], &x_up[0],
          function_->dim, 25000,
          rng_, veg_state,
          &result, &abserr );
        grid_.grid_prepared = true;
      }
      Information( "Finished the Vegas warm-up." );
      //----- integration
      unsigned short it_chisq = 0;
      do {
        res = gsl_monte_vegas_integrate( function_.get(),
          &x_low[0], &x_up[0],
          function_->dim, 0.2 * input_params_->integrator.ncvg,
          rng_, veg_state,
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
        rng_, mis_state,
        &result, &abserr );
    }

    //--- clean integrator
    if      ( algorithm == Plain ) gsl_monte_plain_free( pln_state );
    else if ( algorithm == Vegas ) gsl_monte_vegas_free( veg_state );
    else if ( algorithm == MISER ) gsl_monte_miser_free( mis_state );

#ifdef TIMING_ANALYSIS
    std::cout << "|steps: ";
    for ( const auto& st : steps )
    std::cout << "|" << st/num_points;
    std::cout << std::endl;
#endif

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
    unsigned long i = 0;
    if ( callback )
      callback( *input_params_->generation.last_event, i );
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

  ThreadArgs::ThreadArgs( gsl_rng* rng, gsl_monte_function* function, GridParameters* grid ) :
    ps_bin_( 0 ), rng_( rng ), function_( function ), grid_( grid )
  {
    if ( function )
      params_ = (Parameters*)function->params;
  }

  void
  ThreadArgs::initialise()
  {
    Information( Form( "Preparing the grid (%d points) for the generation of unweighted events.", params_->integrator.npoints ) );
    // Variables for debugging
    std::ostringstream os;
    if ( Logger::get().level >= Logger::Debug )
      Debugging( Form( "Maximum weight = %d", params_->generation.maxgen ) );

    const unsigned int max = pow( mbin_, function_->dim );
    const double inv_npoin = 1./params_->integrator.npoints;

    if ( function_->dim > max_dimensions_ )
      FatalError( Form( "Number of dimensions to integrate exceeds the maximum number, %d", max_dimensions_ ) );

    grid_->nm = std::vector<int>( max, 0 );
    grid_->f_max = std::vector<double>( max, 0. );
    grid_->n = std::vector<int>( function_->dim, 0 );

    std::vector<double> x( function_->dim, 0. );

    params_->generation.ngen = 0;
    params_->setStorage( false );

    // ...
    double sum = 0., sum2 = 0., sum2p = 0.;

    //--- main loop
    for ( unsigned int i = 0; i < max; ++i ) {
      int jj = i;
      for ( unsigned int j = 0; j < function_->dim; ++j ) {
        int jjj = jj*inv_mbin_;
        grid_->n[j] = jj-jjj*mbin_;
        jj = jjj;
      }
      double fsum = 0., fsum2 = 0.;
      for ( unsigned int j = 0; j < params_->integrator.npoints; ++j ) {
        for ( unsigned int k = 0; k < function_->dim; ++k )
          x[k] = ( uniform() + grid_->n[k] ) * inv_mbin_;
        const double z = eval( x );
        grid_->f_max[i] = std::max( grid_->f_max[i], z );
        fsum += z;
        fsum2 += z*z;
      }
      const double av = fsum*inv_npoin, av2 = fsum2*inv_npoin, sig2 = av2 - av*av;
      sum += av;
      sum2 += av2;
      sum2p += sig2;
      grid_->f_max_global = std::max( grid_->f_max_global, grid_->f_max[i] );

      if ( Logger::get().level >= Logger::DebugInsideLoop ) {
        const double sig = sqrt( sig2 );
        const double eff = ( grid_->f_max[i] != 0. ) ? grid_->f_max[i]/av : 1.e4;
        os.str(""); for ( unsigned int j = 0; j < function_->dim; ++j ) { os << grid_->n[j]; if ( j != function_->dim-1 ) os << ", "; }
        DebuggingInsideLoop( Form( "In iteration #%d:\n\t"
                                   "av   = %f\n\t"
                                   "sig  = %f\n\t"
                                   "fmax = %f\n\t"
                                   "eff  = %f\n\t"
                                   "n = (%s)",
                                   i, av, sig, grid_->f_max[i], eff, os.str().c_str() ) );
      }
    } // end of main loop

    sum = sum/max;
    sum2 = sum2/max;
    sum2p = sum2p/max;

    if ( Logger::get().level >= Logger::Debug ) {
      const double sig = sqrt( sum2-sum*sum ), sigp = sqrt( sum2p );

      double eff1 = 0.;
      for ( unsigned int i = 0; i < max; ++i )
        eff1 += ( grid_->f_max[i] / ( max*sum ) );
      const double eff2 = grid_->f_max_global/sum;

      Debugging( Form( "Average function value     =  sum   = %f\n\t"
                       "Average function value**2  =  sum2  = %f\n\t"
                       "Overall standard deviation =  sig   = %f\n\t"
                       "Average standard deviation =  sigp  = %f\n\t"
                       "Maximum function value     = ffmax  = %f\n\t"
                       "Average inefficiency       =  eff1  = %f\n\t"
                       "Overall inefficiency       =  eff2  = %f\n\t",
                       sum, sum2, sig, sigp, grid_->f_max_global, eff1, eff2 ) );
    }
    grid_->gen_prepared = true;
    Information( "Grid prepared! Now launching the production." );
  }

  bool
  ThreadArgs::generate()
  {
    if ( !grid_->gen_prepared )
      initialise();

    const unsigned int max = pow( mbin_, function_->dim );

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
        int jjj = jj / mbin_;
        grid_->n[i] = jj - jjj * mbin_;
        x[i] = ( uniform() + grid_->n[i] ) / mbin_;
        jj = jjj;
      }

      // Get weight for selected x value
      weight = eval( x );
      if ( weight <= 0. ) continue;
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
  ThreadArgs::correctionCycle( std::vector<double>& x, bool& has_correction )
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
        xtmp.emplace_back( ( uniform() + grid_->n[k] ) * inv_mbin_ );
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
  ThreadArgs::eval( const std::vector<double>& x )
  {
    return function_->f( (double*)&x[0], function_->dim, (void*)params_ );
  }

  bool
  ThreadArgs::storeEvent( const std::vector<double>& x )
  {
    params_->setStorage( true );
    const double weight = eval( x );
    params_->setStorage( false );

    if ( weight <= 0. )
      return false;

    params_->generation.ngen += 1;
    if ( params_->generation.ngen % params_->generation.gen_print_every == 0 ) {
      Information( Form( "Generated events: %d", params_->generation.ngen ) );
      params_->generation.last_event->dump();
    }
    return true;
  }
}

