#include "Vegas.h"

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
  
  Debug(Form("Number of integration dimensions: %d\n\t"
             "Number of iterations:             %d\n\t"
             "Number of function calls:         %d", dim_, inParam_->itvg, inParam_->ncvg));

  fFunction = new gsl_monte_function;
  fFunction->f = f_;
  fFunction->dim = dim_;
  fFunction->params = (void*)(inParam_);
  fNumConverg = inParam_->ncvg;
  fNumIter = inParam_->itvg; 
}

Vegas::~Vegas()
{
  Debug("Destructor called");
  
  delete[] fXlow;
  delete[] fXup;
  delete[] fNm;
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
  gsl_monte_vegas_state* state;
  gsl_rng_env_setup();
  rng_type = gsl_rng_default;
  int veg_res;
  
  double res, err;

  // Prepare Vegas
  rng = gsl_rng_alloc(rng_type);
  state = gsl_monte_vegas_alloc(fFunction->dim);
  
  // Launch Vegas
  /// Warmup (prepare the grid)
  if (!fGridPrepared) {
    veg_res = gsl_monte_vegas_integrate(fFunction, fXlow, fXup, fFunction->dim, 10000, rng, state, &res, &err);
    fGridPrepared = true;
  }
  /// Integration
  for (unsigned int i=0; i<fNumIter; i++) {
    veg_res = gsl_monte_vegas_integrate(fFunction, fXlow, fXup, fFunction->dim, fNumConverg/5, rng, state, &res, &err);
    std::cout << Form(">> Iteration %2d: average = %8.4f   sigma = %8.4f   chi2 = %4.3f", i+1, res, err, gsl_monte_vegas_chisq(state)) << std::endl;
  }
  
  // Clean Vegas
  gsl_monte_vegas_free(state);
  gsl_rng_free(rng);
  
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

  i = 0;

  Info(Form("%d events will be generated", fInputParameters->maxgen));
  while (i<fInputParameters->maxgen) {
    if (this->GenerateOneEvent()) i++;
  }
  Info(Form("%d events generated", i));
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
      fNm[fJ] += 1;
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
    fCorrec = (fNm[fJ]-1.)*fmax_diff/fFGlobalMax-1.;
  }
  else {
    fmax_old = fFmax[fJ];
    fFmax[fJ] = weight;
    fmax_diff = weight-fmax_old;
    fFGlobalMax = weight;
    fCorrec = (fNm[fJ]-1.)*fmax_diff/fFGlobalMax*weight/fFGlobalMax-1.;
  }
  
  Debug(Form("Correc.: %f, j = %d", fCorrec, fJ));
  
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
  
  Debug(Form("Correction cycles are started.\n\t"
                  "j = %f"
                  "correc = %f"
                  "corre2 = %f", fJ, fCorrec2));
  
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
      fCorrec = (fNm[fJ]-1.)*fmax_diff/fFGlobalMax-fCorrec2;
    }
    else {
      fFGlobalMax = fmax2;
      fCorrec = (fNm[fJ]-1.)*fmax_diff/fFGlobalMax*fmax2/fFGlobalMax-fCorrec2;
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
  
  if (fInputParameters->ngen%1000==0) {
    Debug(Form("Generated events: %d", fInputParameters->ngen));
  }
  
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

  // Variables for debugging
  double eff, eff1, eff2;
  double sig, sigp;
  std::ostringstream os;
  if (Logger::GetInstance()->Level>=Logger::Debug) {
    Debug(Form("MaxGen = %d", fInputParameters->maxgen));
  }

  fNm = new int[20000];
  fFmax = new double[20000];
  fN = new int[fFunction->dim];

  fInputParameters->ngen = 0;

  // ...
  sum = 0.;
  sum2 = 0.;
  sum2p = 0.;
  max = pow(fMbin, fFunction->dim);

  for (int i=0; i<max; i++) {
    fNm[i] = 0;
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

    if (Logger::GetInstance()->Level>=Logger::Debug) {
      sig = sqrt(sig2);
      eff = 1.e4;
      if (fFmax[i]!=0.) eff = fFmax[i]/av;
      os.str("");
      for (unsigned int j=0; j<fFunction->dim; j++) { os << n[j]; if (j!=fFunction->dim-1) os << ", "; }
      DebugInsideLoop(Form("In iteration #%d:\n\t"
	                              "av   = %f\n\t"
	                              "sig  = %f\n\t"
                                "fmax = %f\n\t"
                                "eff  = %f\n\t"
                                "n = (%s)",
                                i, av, sig, fFmax[i], eff, os.str().c_str()));
    }
  }

  sum = sum/max;
  sum2 = sum2/max;
  sum2p = sum2p/max;

  if (Logger::GetInstance()->Level>=Logger::Debug) {
    sig = sqrt(sum2-pow(sum, 2));
    sigp = sqrt(sum2p);
    
    eff1 = 0.;
    for (int i=0; i<max; i++) eff1 += (fFmax[i]/(max*sum));
    eff2 = fFGlobalMax/sum;
    
    Debug(Form("Average function value     =  sum   = %f\n\t"
                    "Average function value     =  sum   = %f\n\t"
                    "Average function value**2  =  sum2  = %f\n\t"
                    "Overall standard deviation =  sig   = %f\n\t"
                    "Average standard deviation =  sigp  = %f\n\t"
                    "Maximum function value     = ffmax  = %f\n\t"
                    "Average inefficiency       =  eff1  = %f\n\t"
                    "Overall inefficiency       =  eff2  = %f\n\t"
                    "eff = %f", sum, sum2, sig, sigp, fFGlobalMax, eff1, eff2, eff));
  }
}

