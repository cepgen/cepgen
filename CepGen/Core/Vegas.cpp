#include "Vegas.h"

using namespace CepGen;

Vegas::Vegas( const unsigned int dim, double f_( double*, size_t, void* ), Parameters* param ) :
  mbin_( 3 ),
  j_( 0 ), correc_( 0. ), correc2_( 0. ),
  input_params_( param ),
  grid_prepared_( false ), gen_prepared_( false ),
  f_max_( 0 ), f_max2_( 0. ), f_max_diff_( 0. ), f_max_old_( 0. ), f_max_global_( 0. ),
  n_( 0 ), nm_( NULL ), function_( 0 ), x_( 0 )
{
  x_low_ = new double[dim];
  x_up_ = new double[dim];
  
  for ( unsigned int i=0; i<dim; i++ ) {
    x_low_[i] = 0.;
    x_up_[i] = 1.;
  }
  
  Debugging( Form( "Number of integration dimensions: %d\n\t"
                   "Number of iterations:             %d\n\t"
                   "Number of function calls:         %d", dim, param->itvg, param->ncvg ) );

  function_ = new gsl_monte_function;
  function_->f = f_;
  function_->dim = dim;
  function_->params = (void*)param;
  num_converg_ = param->ncvg;
  num_iter_ = param->itvg;
std::cout << __PRETTY_FUNCTION__ << ": " << function_->dim << std::endl;
}

Vegas::~Vegas()
{
  if ( x_low_ ) delete[] x_low_;
  if ( x_up_ ) delete[] x_up_;
  if ( f_max_ ) delete[] f_max_;
  if ( nm_ ) delete[] nm_;
  if ( n_ ) delete[] n_;
  if ( function_ ) delete function_;
  if ( x_ ) delete[] x_;
std::cout << __PRETTY_FUNCTION__ << ": " << function_->dim << std::endl;
}

int
Vegas::integrate( double& result, double& abserr )
{
  // Initialise the random number generator
  const gsl_rng_type* rng_type;
  gsl_rng* rng;
  gsl_monte_vegas_state* state;
  gsl_rng_env_setup();
  rng_type = gsl_rng_default;
  int veg_res;
  
  // Prepare the integration coordinates
  if ( x_ ) delete x_;
  x_ = new double[function_->dim];

  // Prepare Vegas
  rng = gsl_rng_alloc( rng_type );
  state = gsl_monte_vegas_alloc( function_->dim );
  
  // Launch Vegas
  /// Warmup (prepare the grid)
  if ( !grid_prepared_ ) {
    veg_res = gsl_monte_vegas_integrate( function_, x_low_, x_up_, function_->dim, 10000, rng, state, &result, &abserr );
    grid_prepared_ = true;
  }
  /// Integration
  for ( unsigned int i=0; i<num_iter_; i++ ) {
    veg_res = gsl_monte_vegas_integrate( function_, x_low_, x_up_, function_->dim, 0.2*num_converg_, rng, state, &result, &abserr );
    Print( Form( ">> Iteration %2d: average = %10.6f   sigma = %10.6f   chi2 = %4.3f", i+1, result, abserr, gsl_monte_vegas_chisq( state ) ) );
  }
  
  // Clean Vegas
  gsl_monte_vegas_free( state );
  gsl_rng_free( rng );
  
  return veg_res;
}

void
Vegas::generate()
{
  std::ofstream of;
  std::string fn;
  
  this->setGen();

  Information( Form( "%d events will be generated", input_params_->maxgen ) );

  unsigned int i = 0;
  while ( i<input_params_->maxgen ) {
    if ( this->generateOneEvent() ) i++;
  }
  Information( Form( "%d events generated", i ) );
}

bool
Vegas::generateOneEvent()
{
  // Inherited from GMUGNA
  int jj, jjj;
  
  if ( !gen_prepared_ ) setGen();

  const unsigned int ndim = function_->dim,
                     max = pow( mbin_, ndim );

  // Correction cycles are started
  if ( j_!=0 ) {
    bool has_correction = false;
    while ( !correctionCycle( x_, has_correction ) ) {}
    if ( has_correction ) return storeEvent( x_ );
  }

  double weight;
  double y = -1.;

  // Normal generation cycle
  // Select a Vegas bin and reject if fmax is too little
  do {
    do {
      // ...
      j_ = (double)rand()/RAND_MAX * max;
      y = (double)rand()/RAND_MAX * f_max_global_;
      nm_[j_] += 1;
    } while ( y>f_max_[j_] );
    // Select x values in this Vegas bin
    jj = j_;
    for (unsigned int i=0; i<ndim; i++) {
      jjj = jj / mbin_;
      n_[i] = jj - jjj * mbin_;
      x_[i] = ( (double)rand()/RAND_MAX + n_[i] ) / mbin_;
      jj = jjj;
    }
    
    // Get weight for selected x value
    weight = F( x_ );
    
    // Eject if weight is too low
    //if (y>weight) {
      //std::cout << "ERROR: y>weight => " << y << ">" << weight << ", " << j_ << std::endl;
      //_force_correction = false;
      //_force_returnto1 = true;
      //return this->generateOneEvent();
      //goto line1;
    //}
  } while ( y>weight );

  if ( weight<=f_max_[j_] ) j_ = 0;
  // Init correction cycle if weight is higher than fmax or ffmax
  else if ( weight<=f_max_global_ ) {
    f_max_old_ = f_max_[j_];
    f_max_[j_] = weight;
    f_max_diff_ = weight-f_max_old_;
    correc_ = ( nm_[j_] - 1. ) * f_max_diff_ / f_max_global_ - 1.;
  }
  else {
    f_max_old_ = f_max_[j_];
    f_max_[j_] = weight;
    f_max_diff_ = weight-f_max_old_;
    f_max_global_ = weight;
    correc_ = ( nm_[j_] - 1. ) * f_max_diff_ / f_max_global_ * weight / f_max_global_ - 1.;
  }
  
  Debugging( Form( "Correc.: %f, j = %d", correc_, j_ ) );
  
  // Return with an accepted event
  if ( weight>0. ) return storeEvent( x_ );
  return false;
}

bool
Vegas::correctionCycle( double* x_, bool& has_correction )
{
  double weight;
  const unsigned int ndim = function_->dim;
  
  Debugging( Form( "Correction cycles are started.\n\t"
                   "j = %f"
                   "correc = %f"
                   "corre2 = %f", j_, correc2_ ) );
  
  if ( correc_>=1. ) correc_ -= 1.;
  if ( (double)rand()/RAND_MAX<correc_ ) {
    correc_ = -1.;
    // Select x values in Vegas bin
    for ( unsigned int k=0; k<ndim; k++ ) {
      x_[k] = ( (double)rand()/RAND_MAX + n_[k] ) / mbin_;
    }
    // Compute weight for x value
    weight = F( x_ );
    // Parameter for correction of correction
    if ( weight>f_max_[j_] ) {
      if ( weight>f_max2_ ) f_max2_ = weight;
      correc2_ -= 1.;
      correc_ += 1.;
    }
    // Accept event
    if ( weight>=f_max_diff_*(double)rand()/RAND_MAX + f_max_old_ ) { // FIXME!!!!
      //Error("Accepting event!!!");
      //return storeEvent(x);
      std::copy( x_, x_+ndim, x_ );
      has_correction = true;
      return true;
    }
    return false;
  }
  // Correction if too big weight is found while correction
  // (All your bases are belong to us...)
  if ( f_max2_>f_max_[j_] ) {
    f_max_old_ = f_max_[j_];
    f_max_[j_] = f_max2_;
    f_max_diff_ = f_max2_-f_max_old_;
    if ( f_max2_<f_max_global_ ) {
      correc_ = ( nm_[j_] - 1. ) * f_max_diff_ / f_max_global_ - correc2_;
    }
    else {
      f_max_global_ = f_max2_;
      correc_ = ( nm_[j_] - 1. ) * f_max_diff_ / f_max_global_ * f_max2_ / f_max_global_ - correc2_;
    }
    correc2_ = 0.;
    f_max2_ = 0.;
    return false;
  }
  return true;
}

bool
Vegas::storeEvent( double *x_ )
{
  input_params_->store = true;
  F( x_ );
  input_params_->ngen += 1;
  input_params_->store = false;
  
  if ( input_params_->ngen%1000==0 ) {
    Debugging( Form( "Generated events: %d", input_params_->ngen ) );
  }
  
  return true;
}

void
Vegas::setGen()
{
  int jj, jjj;
  double sum, sum2, sum2p;
  int n[15];
  double fsum, fsum2;
  double z;
  double sig2;
  double av, av2;

  // Variables for debugging
  double eff, eff1, eff2;
  double sig, sigp;
  std::ostringstream os;
  if ( Logger::GetInstance()->Level>=Logger::Debug ) {
    Debugging( Form( "MaxGen = %d", input_params_->maxgen ) );
  }
  
  const unsigned int ndim = function_->dim,
                     max = pow( mbin_, ndim ),
                     npoin = input_params_->npoints;

  nm_ = new int[max];
  f_max_ = new double[max];
  n_ = new int[ndim];

  input_params_->ngen = 0;

  // ...
  sum = 0.;
  sum2 = 0.;
  sum2p = 0.;

  for ( unsigned int i=0; i<max; i++ ) {
    nm_[i] = 0;
    f_max_[i] = 0.;
  }

  for ( unsigned int i=0; i<max; i++ ) {
    jj = i;
    for ( unsigned int j=0; j<ndim; j++ ) {
      jjj = jj/mbin_;
      n[j] = jj-jjj*mbin_;
      jj = jjj;
    }
    fsum = fsum2 = 0.;
    for ( unsigned int j=0; j<npoin; j++ ) {
      for ( unsigned int k=0; k<ndim; k++ ) {
        x_[k] = ( (double)rand()/RAND_MAX + n[k] ) / mbin_;
      }
      z = F( x_ );
      if ( z>f_max_[i] ) f_max_[i] = z;
      fsum += z;
      fsum2 += z*z;
    }
    av = fsum/npoin;
    av2 = fsum2/npoin;
    sig2 = av2 - av*av;
    sum += av;
    sum2 += av2;
    sum2p += sig2;
    if ( f_max_[i]>f_max_global_ ) f_max_global_ = f_max_[i];

    if ( Logger::GetInstance()->Level>=Logger::Debug ) {
      sig = sqrt( sig2 );
      eff = 1.e4;
      if ( f_max_[i]!=0. ) eff = f_max_[i]/av;
      os.str(""); for ( unsigned int j=0; j<ndim; j++ ) { os << n[j]; if ( j!=ndim-1 ) os << ", "; }
      DebuggingInsideLoop( Form( "In iteration #%d:\n\t"
                                 "av   = %f\n\t"
                                 "sig  = %f\n\t"
                                 "fmax = %f\n\t"
                                 "eff  = %f\n\t"
                                 "n = (%s)",
                                 i, av, sig, f_max_[i], eff, os.str().c_str() ) );
    }
  }

  sum = sum/max;
  sum2 = sum2/max;
  sum2p = sum2p/max;

  if ( Logger::GetInstance()->Level>=Logger::Debug ) {
    sig = sqrt( sum2-sum*sum );
    sigp = sqrt( sum2p );
    
    eff1 = 0.;
    for ( unsigned int i=0; i<max; i++ ) eff1 += ( f_max_[i] / ( max*sum ) );
    eff2 = f_max_global_/sum;
    
    Debugging( Form( "Average function value     =  sum   = %f\n\t"
                     "Average function value     =  sum   = %f\n\t"
                     "Average function value**2  =  sum2  = %f\n\t"
                     "Overall standard deviation =  sig   = %f\n\t"
                     "Average standard deviation =  sigp  = %f\n\t"
                     "Maximum function value     = ffmax  = %f\n\t"
                     "Average inefficiency       =  eff1  = %f\n\t"
                     "Overall inefficiency       =  eff2  = %f\n\t"
                     "eff = %f",
                     sum, sum2, sig, sigp, f_max_global_, eff1, eff2, eff ) );
  }
  gen_prepared_ = true;
}

