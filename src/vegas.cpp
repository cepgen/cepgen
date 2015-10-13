#include "vegas.h"

Vegas::Vegas(const int dim_, double f_(double*,size_t,void*), Parameters* inParam_) :
  fMbin(3),
  fJ(0), fCorrec(0.), fCorrec2(0.),
  fInputParameters(inParam_),
  fGridPrepared(false), fGenerationPrepared(false),
  fFmax(0), fFGlobalMax(0.), fN(0)
{
  fXlow = new double[dim_];
  fXup = new double[dim_];
  
  for (int i=0; i<dim_; i++) {
    fXlow[i] = 0.;
    fXup[i] = 1.;
  }
  
#ifdef DEBUG
  std::cout << __PRETTY_FUNCTION__ << " [DEBUG]"
            << "\n  Number of integration dimensions : " << dim_
            << "\n  Number of iterations : " << inParam_->itvg
            << "\n  Number of function calls : " << inParam->ncvg
            << std::endl;
#endif

  fFunction = new gsl_monte_function;
  fFunction->f = f_;
  fFunction->dim = dim_;
  fFunction->params = (void*)(inParam_);
  fNumConverg = inParam_->ncvg;
  fNumIter = inParam_->itvg; 
}

Vegas::~Vegas()
{
#ifdef DEBUG
  std::cout << __PRETTY_FUNCTION__ << " [DEBUG] Destructor called" << std::endl;
#endif
  delete[] fXlow;
  delete[] fXup;
  delete[] _nm;
  if (fFmax) delete[] fFmax;
  if (fN) delete[] fN;
  delete fFunction;
}

int
Vegas::Integrate(double *result_, double *abserr_)
{
  // Initialise the random number generator
  const gsl_rng_type* rng_type;
  gsl_rng* rng;
  gsl_rng_env_setup();
  rng_type = gsl_rng_default;
  rng = gsl_rng_alloc(rng_type);
  int veg_res;
  
  double res, err;
  
  // Launch Vegas
  gsl_monte_vegas_state* state = gsl_monte_vegas_alloc(fFunction->dim);
  // Vegas warmup
  if (!fGridPrepared) {
    veg_res = gsl_monte_vegas_integrate(fFunction, fXlow, fXup, fFunction->dim, 10000, rng, state, &res, &err);
    fGridPrepared = true;
  }
  for (unsigned int i=0; i<fNumIter; i++) {
    veg_res = gsl_monte_vegas_integrate(fFunction, fXlow, fXup, fFunction->dim, fNumConverg/5, rng, state, &res, &err);
    std::cout << "--> iteration " 
      	      << std::setfill(' ') << std::setw(2) << (i+1) << " : "
	            << "average = " << std::setprecision(5) << std::setw(14) << res
	            << "sigma = " << std::setprecision(5) << std::setw(14) << err
      	      << "chi2 = " << gsl_monte_vegas_chisq(state)
              << std::endl;
  }
  
  *result_ = res;
  *abserr_ = err;
  return veg_res;
}

void
Vegas::Generate()
{
  std::ofstream of;
  std::string fn;
  int i;
  
  this->SetGen();
  std::cout << __PRETTY_FUNCTION__ << " " << fInputParameters->maxgen << " events will be generated" << std::endl;
  i = 0;
  while (i<fInputParameters->maxgen) {
    if (this->GenerateOneEvent()) i++;
  }
  std::cout << __PRETTY_FUNCTION__ << " " << i << " events generated" << std::endl;
}

bool
Vegas::GenerateOneEvent()
{
  // Inherited from GMUGNA
  double weight;
  double max;
  double y;
  int jj, jjj;
  double x[fFunction->dim];
  double fmax_old, fmax_diff;
  
  if (!fGenerationPrepared) {
    this->SetGen();
    fGenerationPrepared = true;
  }

  y = -1.;
  max = pow(fMbin, fFunction->dim);

  // Correction cycles are started
  if (fJ!=0) {
    while (CorrectionCycle()) {;}
  }

  // Normal generation cycle
  // Select a Vegas bin and reject if fmax is too little
  //line1:
  //double* Vegas::SelectBin() //FIXME need to implement it this way instead of these bloody goto...!
  do {
    do {
      // ...
      fJ = (double)rand()/RAND_MAX*max;
      y = (double)rand()/RAND_MAX*fFGlobalMax;
      _nm[fJ] += 1;
    } while (y>fFmax[fJ]);
    // Select x values in this Vegas bin
    jj = fJ;
    for (unsigned int i=0; i<fFunction->dim; i++) {
      jjj = jj/fMbin;
      fN[i] = jj-jjj*fMbin;
      x[i] = ((double)rand()/RAND_MAX+fN[i])/fMbin;
      jj = jjj;
    }
    
    // Get weight for selected x value
    weight = F(x);
    
    // Eject if weight is too low
    //if (y>weight) {
      //std::cout << "ERROR : y>weight => " << y << ">" << weight << ", " << fJ << std::endl;
      //_force_correction = false;
      //_force_returnto1 = true;
      //return this->GenerateOneEvent();
      //goto line1;
    //}
  } while (y>weight);

  if (weight<=fFmax[fJ]) fJ = 0;
  // Init correction cycle if weight is higher than fmax or ffmax
  else if (weight<=fFGlobalMax) {
    fmax_old = fFmax[fJ];
    fFmax[fJ] = weight;
    fmax_diff = weight-fmax_old;
    fCorrec = (_nm[fJ]-1.)*fmax_diff/fFGlobalMax-1.;
  }
  else {
    fmax_old = fFmax[fJ];
    fFmax[fJ] = weight;
    fmax_diff = weight-fmax_old;
    fFGlobalMax = weight;
    fCorrec = (_nm[fJ]-1.)*fmax_diff/fFGlobalMax*weight/fFGlobalMax-1.;
  }
#ifdef DEBUG
  std::cout << __PRETTY_FUNCTION__ << " [DEBUG] correc = " << fCorrec << ", j = " << fJ << std::endl;
#endif
  // Return with an accepted event
  if (weight>0.) return this->StoreEvent(x);
  return false;
}

bool
Vegas::CorrectionCycle()
{
  double x[fFunction->dim];
  double weight;
  double fmax_old = 0., fmax_diff = 0., fmax2 = 0.;
  
#ifdef DEBUG
  std::cout << __PRETTY_FUNCTION__ << " [DEBUG] Correction cycles are started."
      << "\n\tj = " << fJ
      << "\n\tcorrec = " << fCorrec
      << "\n\tcorre2 = " << fCorrec2
      << std::endl;
#endif
  if (fCorrec>=1.) {
    fCorrec -= 1.;
  }
  if ((double)rand()/RAND_MAX<fCorrec) {
    fCorrec = -1.;
    // Select x values in Vegas bin
    for (unsigned int k=0; k<fFunction->dim; k++) {
      x[k] = ((double)rand()/RAND_MAX+fN[k])/fMbin;
    }
    // Compute weight for x value
    weight = F(x);
    // Parameter for correction of correction
    if (weight>fFmax[fJ]) {
      if (weight>fmax2) fmax2 = weight;
      fCorrec2 -= 1.;
      fCorrec += 1.;
    }
    // Accept event
    if (weight>=fmax_diff*(double)rand()/RAND_MAX+fmax_old) { // FIXME!!!!
      return this->StoreEvent(x);
    }
    return false;
  }
  // Correction if too big weight is found while correction
  // (All your bases are belong to us...)
  if (fmax2>fFmax[fJ]) {
    fmax_old = fFmax[fJ];
    fFmax[fJ] = fmax2;
    fmax_diff = fmax2-fmax_old;
    if (fmax2<fFGlobalMax) {
      fCorrec = (_nm[fJ]-1.)*fmax_diff/fFGlobalMax-fCorrec2;
    }
    else {
      fFGlobalMax = fmax2;
      fCorrec = (_nm[fJ]-1.)*fmax_diff/fFGlobalMax*fmax2/fFGlobalMax-fCorrec2;
    }
    fCorrec2 = 0.;
    fmax2 = 0.;
  }
  return false;
}

bool
Vegas::StoreEvent(double *x_)
{
  fInputParameters->store = true;
  F(x_);
  fInputParameters->ngen += 1;
  fInputParameters->store = false;
#ifdef DEBUG
  if (fInputParameters->ngen%1000==0) {
    std::cout << __PRETTY_FUNCTION__ << " Generated events : " << fInputParameters->ngen << std::endl;
  }
#endif
  return true;
}

void
Vegas::SetGen()
{
  int max;
  int jj, jjj;
  double sum, sum2, sum2p;
  int n[10];
  int npoin = fInputParameters->npoints;
  double fsum, fsum2;
  double z;
  double x[fFunction->dim];
  double sig2;
  double av, av2;

//#define DEBUG

#ifdef DEBUG
  double eff, eff1, eff2;
  double sig, sigp;
  
  std::cout << __PRETTY_FUNCTION__ << " [DEBUG] maxgen = " << fInputParameters->maxgen << std::endl;
  fInputParameters->Dump();
#endif

  _nm = new int[20000];
  fFmax = new double[20000];
  fN = new int[fFunction->dim];

  fInputParameters->ngen = 0;

  // ...
  sum = 0.;
  sum2 = 0.;
  sum2p = 0.;
  max = pow(fMbin, fFunction->dim);

  for (int i=0; i<max; i++) {
    _nm[i] = 0;
    fFmax[i] = 0.;
  }

  for (int i=0; i<max; i++) {
    jj = i;
    for (unsigned int j=0; j<fFunction->dim; j++) {
      jjj = jj/fMbin;
      n[j] = jj-jjj*fMbin;
      jj = jjj;
    }
    fsum = fsum2 = 0.;
    for (int j=0; j<npoin; j++) {
      for (unsigned int k=0; k<fFunction->dim; k++) {
        x[k] = ((double)rand()/RAND_MAX+n[k])/fMbin;
      }
      z = F(x);
      if (z>fFmax[i]) fFmax[i] = z;
      fsum += z;
      fsum2 += std::pow(z, 2);
    }
    av = fsum/npoin;
    av2 = fsum2/npoin;
    sig2 = av2-pow(av, 2);
    sum += av;
    sum2 += av2;
    sum2p += sig2;
    if (fFmax[i]>fFGlobalMax) fFGlobalMax = fFmax[i];
#ifdef DEBUG
    sig = sqrt(sig2);
    eff = 1.e4;
    if (fFmax[i]!=0.) eff = fFmax[i]/av;
    //#ifdef DEBUG
    std::cout << __PRETTY_FUNCTION__ << " [DEBUG] in iteration #" << i << " :"
	      << "\n\tav   = " << av
	      << "\n\tsig  = " << sig
	      << "\n\tfmax = " << fFmax[i]
	      << "\n\teff  = " << eff
	      << "\n\tn = (";
    for (unsigned int j=0; j<fFunction->dim; j++) {
      std::cout << n[j];
      if (j!=fFunction->dim-1) std::cout << ", ";
    }
    std::cout << ")" << std::endl;
#endif
  }

  sum = sum/max;
  sum2 = sum2/max;
  sum2p = sum2p/max;

#ifdef DEBUG
  sig = sqrt(sum2-pow(sum, 2));
  sigp = sqrt(sum2p);
  eff1 = 0.;
  for (int i=0; i<max; i++) {
    eff1 += fFmax[i];
  }
  eff1 = eff1/(max*sum);
  eff2 = fFGlobalMax/sum;
  std::cout << __PRETTY_FUNCTION__ << " [DEBUG]"
            << "\n\tAverage function value     =  sum   = " << sum
            << "\n\tAverage function value**2  =  sum2  = " << sum2
            << "\n\tOverall standard deviation =  sig   = " << sig
            << "\n\tAverage standard deviation =  sigp  = " << sigp
            << "\n\tMaximum function value     = ffmax  = " << fFGlobalMax
            << "\n\tAverage inefficiency       =  eff1  = " << eff1 
            << "\n\tOverall inefficiency       =  eff2  = " << eff2 
            << "\n\teff = " << eff 
            << std::endl;
#endif
//#undef DEBUG
}

