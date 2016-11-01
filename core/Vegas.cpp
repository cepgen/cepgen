#include "Vegas.h"

Vegas::Vegas( const unsigned int dim_, double f_( double*, size_t, void* ), Parameters* inParam_ ) :
  fMbin( 3 ),
  fJ( 0 ), fCorrec( 0. ), fCorrec2( 0. ),
  fInputParameters( inParam_ ),
  fGridPrepared( false ), fGenerationPrepared( false ),
  fFmax( 0 ), fFmax2( 0. ), fFmaxDiff( 0. ), fFmaxOld( 0. ), fFGlobalMax( 0. ),
  fN( 0 ), fNm( NULL ), fFunction( 0 ), fX( 0 )
{
  fXlow = new double[dim_];
  fXup = new double[dim_];
  
  for ( unsigned int i=0; i<dim_; i++ ) {
    fXlow[i] = 0.;
    fXup[i] = 1.;
  }
  
  Debugging( Form( "Number of integration dimensions: %d\n\t"
                   "Number of iterations:             %d\n\t"
                   "Number of function calls:         %d", dim_, inParam_->itvg, inParam_->ncvg ) );

  fFunction = new gsl_monte_function;
  fFunction->f = f_;
  fFunction->dim = dim_;
  fFunction->params = (void*)inParam_;
  fNumConverg = inParam_->ncvg;
  fNumIter = inParam_->itvg; 
}

Vegas::~Vegas()
{
  if ( fXlow ) delete[] fXlow;
  if ( fXup ) delete[] fXup;
  if ( fFmax ) delete[] fFmax;
  if ( fNm ) delete[] fNm;
  if ( fN ) delete[] fN;
  if ( fFunction ) delete fFunction;
  if ( fX ) delete[] fX;
}

int
Vegas::Integrate( double *result_, double *abserr_ )
{
  // Initialise the random number generator
  const gsl_rng_type* rng_type;
  gsl_rng* rng;
  gsl_monte_vegas_state* state;
  gsl_rng_env_setup();
  rng_type = gsl_rng_default;
  int veg_res;
  
  double res, err;

  // Prepare the integration coordinates
  if ( fX ) delete fX;
  fX = new double[fFunction->dim];

  // Prepare Vegas
  rng = gsl_rng_alloc( rng_type );
  state = gsl_monte_vegas_alloc( fFunction->dim );
  
  // Launch Vegas
  /// Warmup (prepare the grid)
  if ( !fGridPrepared ) {
    veg_res = gsl_monte_vegas_integrate( fFunction, fXlow, fXup, fFunction->dim, 10000, rng, state, &res, &err );
    fGridPrepared = true;
  }
  /// Integration
  for ( unsigned int i=0; i<fNumIter; i++ ) {
    veg_res = gsl_monte_vegas_integrate( fFunction, fXlow, fXup, fFunction->dim, fNumConverg/5, rng, state, &res, &err );
    std::cout << Form( ">> Iteration %2d: average = %8.4f   sigma = %8.4f   chi2 = %4.3f", i+1, res, err, gsl_monte_vegas_chisq( state ) ) << std::endl;
  }
  
  // Clean Vegas
  gsl_monte_vegas_free( state );
  gsl_rng_free( rng );
  
  *result_ = res;
  *abserr_ = err;
  
  return veg_res;
}

void
Vegas::Generate()
{
  std::ofstream of;
  std::string fn;
  
  this->SetGen();

  Information( Form( "%d events will be generated", fInputParameters->maxgen ) );

  unsigned int i = 0;
  while ( i<fInputParameters->maxgen ) {
    if ( this->GenerateOneEvent() ) i++;
  }
  Information( Form( "%d events generated", i ) );
}

bool
Vegas::GenerateOneEvent()
{
  // Inherited from GMUGNA
  int jj, jjj;
  
  if ( !fGenerationPrepared ) SetGen();

  const unsigned int ndim = fFunction->dim,
                     max = pow( fMbin, ndim );

  // Correction cycles are started
  if ( fJ!=0 ) {
    fHasCorrection = false;
    while ( !CorrectionCycle(fX) ) {;}
    if ( fHasCorrection ) return StoreEvent( fX );
  }

  double weight;
  double y = -1.;

  // Normal generation cycle
  // Select a Vegas bin and reject if fmax is too little
  do {
    do {
      // ...
      fJ = (double)rand()/RAND_MAX * max;
      y = (double)rand()/RAND_MAX * fFGlobalMax;
      fNm[fJ] += 1;
    } while ( y>fFmax[fJ] );
    // Select x values in this Vegas bin
    jj = fJ;
    for (unsigned int i=0; i<ndim; i++) {
      jjj = jj / fMbin;
      fN[i] = jj - jjj * fMbin;
      fX[i] = ( (double)rand()/RAND_MAX + fN[i] ) / fMbin;
      jj = jjj;
    }
    
    // Get weight for selected x value
    weight = F( fX );
    
    // Eject if weight is too low
    //if (y>weight) {
      //std::cout << "ERROR : y>weight => " << y << ">" << weight << ", " << fJ << std::endl;
      //_force_correction = false;
      //_force_returnto1 = true;
      //return this->GenerateOneEvent();
      //goto line1;
    //}
  } while ( y>weight );

  if ( weight<=fFmax[fJ] ) fJ = 0;
  // Init correction cycle if weight is higher than fmax or ffmax
  else if ( weight<=fFGlobalMax ) {
    fFmaxOld = fFmax[fJ];
    fFmax[fJ] = weight;
    fFmaxDiff = weight-fFmaxOld;
    fCorrec = ( fNm[fJ] - 1. ) * fFmaxDiff / fFGlobalMax - 1.;
  }
  else {
    fFmaxOld = fFmax[fJ];
    fFmax[fJ] = weight;
    fFmaxDiff = weight-fFmaxOld;
    fFGlobalMax = weight;
    fCorrec = ( fNm[fJ] - 1. ) * fFmaxDiff / fFGlobalMax * weight / fFGlobalMax - 1.;
  }
  
  Debugging( Form( "Correc.: %f, j = %d", fCorrec, fJ ) );
  
  // Return with an accepted event
  if ( weight>0. ) return StoreEvent( fX );
  return false;
}

bool
Vegas::CorrectionCycle( double* x_ )
{
  double weight;
  const unsigned int ndim = fFunction->dim;
  
  Debugging( Form( "Correction cycles are started.\n\t"
                   "j = %f"
                   "correc = %f"
                   "corre2 = %f", fJ, fCorrec2 ) );
  
  if ( fCorrec>=1. ) fCorrec -= 1.;
  if ( (double)rand()/RAND_MAX<fCorrec ) {
    fCorrec = -1.;
    // Select x values in Vegas bin
    for ( unsigned int k=0; k<ndim; k++ ) {
      fX[k] = ( (double)rand()/RAND_MAX + fN[k] ) / fMbin;
    }
    // Compute weight for x value
    weight = F( fX );
    // Parameter for correction of correction
    if ( weight>fFmax[fJ] ) {
      if ( weight>fFmax2 ) fFmax2 = weight;
      fCorrec2 -= 1.;
      fCorrec += 1.;
    }
    // Accept event
    if ( weight>=fFmaxDiff*(double)rand()/RAND_MAX + fFmaxOld ) { // FIXME!!!!
      //Error("Accepting event!!!");
      //return StoreEvent(x);
      std::copy( fX, fX+ndim, x_ );
      fHasCorrection = true;
      return true;
    }
    return false;
  }
  // Correction if too big weight is found while correction
  // (All your bases are belong to us...)
  if ( fFmax2>fFmax[fJ] ) {
    fFmaxOld = fFmax[fJ];
    fFmax[fJ] = fFmax2;
    fFmaxDiff = fFmax2-fFmaxOld;
    if ( fFmax2<fFGlobalMax ) {
      fCorrec = ( fNm[fJ] - 1. ) * fFmaxDiff / fFGlobalMax - fCorrec2;
    }
    else {
      fFGlobalMax = fFmax2;
      fCorrec = ( fNm[fJ] - 1. ) * fFmaxDiff / fFGlobalMax * fFmax2 / fFGlobalMax - fCorrec2;
    }
    fCorrec2 = 0.;
    fFmax2 = 0.;
    return false;
  }
  return true;
}

bool
Vegas::StoreEvent( double *x_ )
{
  fInputParameters->store = true;
  F( x_ );
  fInputParameters->ngen += 1;
  fInputParameters->store = false;
  
  if ( fInputParameters->ngen%1000==0 ) {
    Debugging( Form( "Generated events: %d", fInputParameters->ngen ) );
  }
  
  return true;
}

void
Vegas::SetGen()
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
    Debugging( Form( "MaxGen = %d", fInputParameters->maxgen ) );
  }
  
  const unsigned int ndim = fFunction->dim,
                     max = pow( fMbin, ndim ),
                     npoin = fInputParameters->npoints;

  fNm = new int[max];
  fFmax = new double[max];
  fN = new int[ndim];

  fInputParameters->ngen = 0;

  // ...
  sum = 0.;
  sum2 = 0.;
  sum2p = 0.;

  for ( unsigned int i=0; i<max; i++ ) {
    fNm[i] = 0;
    fFmax[i] = 0.;
  }

  for ( unsigned int i=0; i<max; i++ ) {
    jj = i;
    for ( unsigned int j=0; j<ndim; j++ ) {
      jjj = jj/fMbin;
      n[j] = jj-jjj*fMbin;
      jj = jjj;
    }
    fsum = fsum2 = 0.;
    for ( unsigned int j=0; j<npoin; j++ ) {
      for ( unsigned int k=0; k<ndim; k++ ) {
        fX[k] = ( (double)rand()/RAND_MAX + n[k] ) / fMbin;
      }
      z = F( fX );
      if ( z>fFmax[i] ) fFmax[i] = z;
      fsum += z;
      fsum2 += z*z;
    }
    av = fsum/npoin;
    av2 = fsum2/npoin;
    sig2 = av2 - av*av;
    sum += av;
    sum2 += av2;
    sum2p += sig2;
    if ( fFmax[i]>fFGlobalMax ) fFGlobalMax = fFmax[i];

    if ( Logger::GetInstance()->Level>=Logger::Debug ) {
      sig = sqrt( sig2 );
      eff = 1.e4;
      if ( fFmax[i]!=0. ) eff = fFmax[i]/av;
      os.str(""); for ( unsigned int j=0; j<ndim; j++ ) { os << n[j]; if ( j!=ndim-1 ) os << ", "; }
      DebuggingInsideLoop( Form( "In iteration #%d:\n\t"
                                 "av   = %f\n\t"
                                 "sig  = %f\n\t"
                                 "fmax = %f\n\t"
                                 "eff  = %f\n\t"
                                 "n = (%s)",
                                 i, av, sig, fFmax[i], eff, os.str().c_str() ) );
    }
  }

  sum = sum/max;
  sum2 = sum2/max;
  sum2p = sum2p/max;

  if ( Logger::GetInstance()->Level>=Logger::Debug ) {
    sig = sqrt( sum2-sum*sum );
    sigp = sqrt( sum2p );
    
    eff1 = 0.;
    for ( unsigned int i=0; i<max; i++ ) eff1 += ( fFmax[i] / ( max*sum ) );
    eff2 = fFGlobalMax/sum;
    
    Debugging( Form( "Average function value     =  sum   = %f\n\t"
                     "Average function value     =  sum   = %f\n\t"
                     "Average function value**2  =  sum2  = %f\n\t"
                     "Overall standard deviation =  sig   = %f\n\t"
                     "Average standard deviation =  sigp  = %f\n\t"
                     "Maximum function value     = ffmax  = %f\n\t"
                     "Average inefficiency       =  eff1  = %f\n\t"
                     "Overall inefficiency       =  eff2  = %f\n\t"
                     "eff = %f",
                     sum, sum2, sig, sigp, fFGlobalMax, eff1, eff2, eff ) );
  }
  fGenerationPrepared = true;
}

