#include "Integrator.h"
#include "CepGen/Parameters.h"
#include "CepGen/Core/Exception.h"

#include <fstream>

namespace CepGen
{
  Integrator::Integrator( const unsigned int dim, double f_( double*, size_t, void* ), Parameters* param ) :
    ps_bin_( 0 ), input_params_( param ), function_( std::unique_ptr<gsl_monte_function>( new gsl_monte_function ) )
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

    Debugging( Form( "Number of integration dimensions: %d\n\t"
                     "Number of iterations:             %d\n\t"
                     "Number of function calls:         %d", dim, input_params_->integrator.itvg, input_params_->integrator.ncvg ) );
  }

  Integrator::~Integrator()
  {
    if ( rng_ ) gsl_rng_free( rng_ );
  }

  int
  Integrator::integrate( double& result, double& abserr )
  {
    int res = -1;
    gsl_monte_vegas_state* veg_state;
    gsl_monte_miser_state* mis_state;
    const Integrator::Type algorithm = input_params_->integrator.type;

    //--- integration bounds
    std::vector<double> x_low( function_->dim, 0. ), x_up( function_->dim, 1. );

    //--- prepare integrator
    if      ( algorithm == Vegas ) veg_state = gsl_monte_vegas_alloc( function_->dim );
    else if ( algorithm == MISER ) mis_state = gsl_monte_miser_alloc( function_->dim );

    if ( algorithm == Vegas ) {
      //----- Vegas warmup (prepare the grid)
      if ( !grid_.grid_prepared ) {
        res = gsl_monte_vegas_integrate( function_.get(), &x_low[0], &x_up[0], function_->dim, 10000, rng_, veg_state, &result, &abserr );
        grid_.grid_prepared = true;
      }
      //----- integration
      for ( unsigned short i = 0; i < input_params_->integrator.itvg; i++ ) {
        res = gsl_monte_vegas_integrate( function_.get(), &x_low[0], &x_up[0], function_->dim, 0.2 * input_params_->integrator.ncvg, rng_, veg_state, &result, &abserr );
        PrintMessage( Form( ">> Iteration %2d: average = %10.6f   sigma = %10.6f   chi2 = %4.3f", i+1, result, abserr, gsl_monte_vegas_chisq( veg_state ) ) );
      }
    }
    //----- integration
    else if ( algorithm == MISER )
      res = gsl_monte_miser_integrate( function_.get(), &x_low[0], &x_up[0], function_->dim, input_params_->integrator.ncvg, rng_, mis_state, &result, &abserr );

    //--- clean integrator
    if      ( algorithm == Vegas ) gsl_monte_vegas_free( veg_state );
    else if ( algorithm == MISER ) gsl_monte_miser_free( mis_state );

    return res;
  }

  void
  Integrator::generate()
  {
    std::ofstream of;
    std::string fn;

    if ( !grid_.gen_prepared ) setGen();

    Information( Form( "%d events will be generated", input_params_->generation.maxgen ) );

    unsigned int i = 0;
    while ( i < input_params_->generation.maxgen ) {
      if ( generateOneEvent() ) i++;
    }
    Information( Form( "%d events generated", i ) );
  }

  bool
  Integrator::generateOneEvent()
  {
    if ( !grid_.gen_prepared ) setGen();

    const unsigned int ndim = function_->dim, max = pow( mbin_, ndim );

    std::vector<double> x( ndim, 0. );

    //--- correction cycles
    
    if ( ps_bin_ != 0 ) {
      bool has_correction = false;
      while ( !correctionCycle( x, has_correction ) ) {}
      if ( has_correction ) return storeEvent( x );
    }

    double weight;
    double y = -1.;

    //--- normal generation cycle

    //----- select a Integrator bin and reject if fmax is too small
    do {
      do {
        // ...
        ps_bin_ = uniform() * max;
        y = uniform() * grid_.f_max_global;
        grid_.nm[ps_bin_] += 1;
      } while ( y > grid_.f_max[ps_bin_] );
      // Select x values in this Integrator bin
      int jj = ps_bin_;
      for ( unsigned int i=0; i<ndim; i++ ) {
        int jjj = jj / mbin_;
        grid_.n[i] = jj - jjj * mbin_;
        x[i] = ( uniform() + grid_.n[i] ) / mbin_;
        jj = jjj;
      }

      // Get weight for selected x value
      weight = F( x );
    } while ( y > weight );

    if ( weight <= grid_.f_max[ps_bin_] ) ps_bin_ = 0;
    // Init correction cycle if weight is higher than fmax or ffmax
    else if ( weight <= grid_.f_max_global ) {
      grid_.f_max_old = grid_.f_max[ps_bin_];
      grid_.f_max[ps_bin_] = weight;
      grid_.f_max_diff = weight-grid_.f_max_old;
      grid_.correc = ( grid_.nm[ps_bin_] - 1. ) * grid_.f_max_diff / grid_.f_max_global - 1.;
    }
    else {
      grid_.f_max_old = grid_.f_max[ps_bin_];
      grid_.f_max[ps_bin_] = weight;
      grid_.f_max_diff = weight-grid_.f_max_old;
      grid_.f_max_global = weight;
      grid_.correc = ( grid_.nm[ps_bin_] - 1. ) * grid_.f_max_diff / grid_.f_max_global * weight / grid_.f_max_global - 1.;
    }

    Debugging( Form( "Correction applied: %f, phase space bin = %d", grid_.correc, ps_bin_ ) );

    // Return with an accepted event
    if ( weight > 0. ) return storeEvent( x );
    return false;
  }

  bool
  Integrator::correctionCycle( std::vector<double>& x, bool& has_correction )
  {
    double weight;
    const unsigned int ndim = function_->dim;

    Debugging( Form( "Correction cycles are started.\n\t"
                     "j = %f"
                     "correc = %f"
                     "corre2 = %f", ps_bin_, grid_.correc2 ) );

    if ( grid_.correc >= 1. ) grid_.correc -= 1.;
    if ( uniform() < grid_.correc ) {
      grid_.correc = -1.;
      std::vector<double> xtmp( ndim, 0. );
      // Select x values in phase space bin
      for ( unsigned int k=0; k<ndim; k++ ) {
        xtmp[k] = ( uniform() + grid_.n[k] ) * inv_mbin_;
      }
      // Compute weight for x value
      weight = F( xtmp );
      // Parameter for correction of correction
      if ( weight > grid_.f_max[ps_bin_] ) {
        if ( weight > grid_.f_max2 ) grid_.f_max2 = weight;
        grid_.correc2 -= 1.;
        grid_.correc += 1.;
      }
      // Accept event
      if ( weight >= grid_.f_max_diff*uniform() + grid_.f_max_old ) { // FIXME!!!!
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
    if ( grid_.f_max2 > grid_.f_max[ps_bin_] ) {
      grid_.f_max_old = grid_.f_max[ps_bin_];
      grid_.f_max[ps_bin_] = grid_.f_max2;
      grid_.f_max_diff = grid_.f_max2-grid_.f_max_old;
      if ( grid_.f_max2 < grid_.f_max_global ) {
        grid_.correc = ( grid_.nm[ps_bin_] - 1. ) * grid_.f_max_diff / grid_.f_max_global - grid_.correc2;
      }
      else {
        grid_.f_max_global = grid_.f_max2;
        grid_.correc = ( grid_.nm[ps_bin_] - 1. ) * grid_.f_max_diff / grid_.f_max_global * grid_.f_max2 / grid_.f_max_global - grid_.correc2;
      }
      grid_.correc2 = 0.;
      grid_.f_max2 = 0.;
      return false;
    }
    return true;
  }

  bool
  Integrator::storeEvent( const std::vector<double>& x )
  {
    input_params_->setStorage( true );
    F( x );
    input_params_->generation.ngen += 1;
    input_params_->setStorage( false );
    if ( input_params_->generation.ngen % input_params_->generation.gen_print_every == 0 ) {
      Debugging( Form( "Generated events: %d", input_params_->generation.ngen ) );
      input_params_->generation.last_event->dump();
    }
    return true;
  }

  void
  Integrator::setGen()
  {
    Information( Form( "Preparing the grid for the generation of unweighted events: %d points", input_params_->integrator.npoints ) );
    // Variables for debugging
    std::ostringstream os;
    if ( Logger::get().level >= Logger::Debug ) {
      Debugging( Form( "MaxGen = %d", input_params_->generation.maxgen ) );
    }

    const unsigned int ndim = function_->dim,
                       max = pow( mbin_, ndim ),
                       npoin = input_params_->integrator.npoints;
    const double inv_npoin = 1./npoin;

    if ( ndim > max_dimensions_ ) {
      FatalError( Form( "Number of dimensions to integrate exceed the maximum number, %d", max_dimensions_ ) );
    }

    grid_.nm = std::vector<int>( max, 0 );
    grid_.f_max = std::vector<double>( max, 0. );
    grid_.n = std::vector<int>( ndim, 0 );

    std::vector<double> x( ndim, 0. );

    input_params_->generation.ngen = 0;

    // ...
    double sum = 0., sum2 = 0., sum2p = 0.;

    //--- main loop
    for ( unsigned int i=0; i<max; i++ ) {
      int jj = i;
      for ( unsigned int j=0; j<ndim; j++ ) {
        int jjj = jj*inv_mbin_;
        grid_.n[j] = jj-jjj*mbin_;
        jj = jjj;
      }
      double fsum = 0., fsum2 = 0.;
      for ( unsigned int j=0; j<npoin; j++ ) {
        for ( unsigned int k=0; k<ndim; k++ ) {
          x[k] = ( uniform() + grid_.n[k] ) * inv_mbin_;
        }
        const double z = F( x );
        grid_.f_max[i] = std::max( grid_.f_max[i], z );
        fsum += z;
        fsum2 += z*z;
      }
      const double av = fsum*inv_npoin, av2 = fsum2*inv_npoin, sig2 = av2 - av*av;
      sum += av;
      sum2 += av2;
      sum2p += sig2;
      grid_.f_max_global = std::max( grid_.f_max_global, grid_.f_max[i] );

      if ( Logger::get().level >= Logger::DebugInsideLoop ) {
        const double sig = sqrt( sig2 );
        const double eff = ( grid_.f_max[i] != 0. ) ? grid_.f_max[i]/av : 1.e4;
        os.str(""); for ( unsigned int j=0; j<ndim; j++ ) { os << grid_.n[j]; if ( j != ndim-1 ) os << ", "; }
        DebuggingInsideLoop( Form( "In iteration #%d:\n\t"
                                   "av   = %f\n\t"
                                   "sig  = %f\n\t"
                                   "fmax = %f\n\t"
                                   "eff  = %f\n\t"
                                   "n = (%s)",
                                   i, av, sig, grid_.f_max[i], eff, os.str().c_str() ) );
      }
    } // end of main loop

    sum = sum/max;
    sum2 = sum2/max;
    sum2p = sum2p/max;

    if ( Logger::get().level >= Logger::Debug ) {
      const double sig = sqrt( sum2-sum*sum ), sigp = sqrt( sum2p );

      double eff1 = 0.;
      for ( unsigned int i=0; i<max; i++ ) eff1 += ( grid_.f_max[i] / ( max*sum ) );
      const double eff2 = grid_.f_max_global/sum;

      Debugging( Form( "Average function value     =  sum   = %f\n\t"
                       "Average function value**2  =  sum2  = %f\n\t"
                       "Overall standard deviation =  sig   = %f\n\t"
                       "Average standard deviation =  sigp  = %f\n\t"
                       "Maximum function value     = ffmax  = %f\n\t"
                       "Average inefficiency       =  eff1  = %f\n\t"
                       "Overall inefficiency       =  eff2  = %f\n\t",
                       sum, sum2, sig, sigp, grid_.f_max_global, eff1, eff2 ) );
    }
    grid_.gen_prepared = true;
    Information( "Grid prepared! Now launching the production." );
  }

  std::ostream&
  operator<<( std::ostream& os, const Integrator::Type& type )
  {
    switch ( type ) {
      case Integrator::Vegas: return os << "Vegas";
      case Integrator::MISER: return os << "MISER";
    }
    return os;
  }
}

