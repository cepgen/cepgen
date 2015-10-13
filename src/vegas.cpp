#include "vegas.h"

Vegas::Vegas(const int dim_, double f_(double*,size_t,void*), Parameters* inParam_) :
  fMbin(3), fStateSamples(0),
  fJ(0), fCorrec(0.), fCorrec2(0.),
  fInputParameters(inParam_),
  fGridPrepared(false), fGenerationPrepared(false),
  fFmax(0), fFGlobalMax(0.), fN(0)
  /*fFunction->dim(dim_), fStateBins(50),
  _nTreatCalls(0), fMbin(3),
  _fmax2(0.), _fmdiff(0.), _fmold(0.),
  fJ(0),
  fStateMode(1), _acc(1.e-4), _alph(1.5), fStateSamples(0)*/
{
  /* x content :
      0 = t1 mapping
      1 = t2 mapping
      2 = s2 mapping
      3 = yy4 definition
      4 = w4 mapping
      5 = xx6 definition
      6 = phicm6 definition
    ( 7 = xq, wx mappings   ) <- single- and double-dissociative only
  */

  fXlow = new double[dim_];
  fXup = new double[dim_];
  /*for (int i=0; i<fMaxNbins; i++) {
    fCoord[i] = new double[dim_];
    fValue[i] = new double[dim_];
    _di[i] = new double[dim_];
  }
  // ...
  // ...*/
  
  for (int i=0; i<dim_; i++) {
    fXlow[i] = 0.;
    fXup[i] = 1.;
  }
  
#ifdef DEBUG
  std::cout << "[Vegas::Vegas] [DEBUG]"
            << "\n  Number of integration dimensions : " << dim_
            << "\n  Number of iterations : " << inParam_->itvg
            << "\n  Number of function calls : " << inParam->ncvg
            << std::endl;
#endif

  /*for (unsigned int j=0; j<fMaxNbins; j++) {
    for (unsigned int i=0; i<fFunction->dim; i++) {
      fValue[j][i] = _di[j][i] = fCoord[j][i] = 0.;
    }
  }*/

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
  std::cout << "[Vegas::~Vegas] [DEBUG] Destructor called" << std::endl;
#endif
  delete[] fXlow;
  delete[] fXup;
  delete[] _nm;
  if (fFmax) delete[] fFmax;
  if (fN) delete[] fN;
  /*for (int i=0; i<fMaxNbins; i++) {
    delete[] fCoord[i];
    delete[] fValue[i];
    delete[] _di[i];
  }*/
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
  
  std::cout << "Will launch vegas on " << fFunction->dim << " dimensions" << std::endl;
  
  // Launch Vegas
  gsl_monte_vegas_state* state = gsl_monte_vegas_alloc(fFunction->dim);
  // Vegas warmup
  if (!fGridPrepared) {
    veg_res = gsl_monte_vegas_integrate(fFunction, fXlow, fXup, fFunction->dim, 10000, rng, state, &res, &err);
    fGridPrepared = true;
  }
  for (unsigned int i=0; i<fNumIter; i++) {
    veg_res = gsl_monte_vegas_integrate(fFunction, fXlow, fXup, fFunction->dim, fNumConverg/5, rng, state, &res, &err);
    printf ("result = % .6f sigma = % .6f chisq/dof = %.1f\n", res, err, gsl_monte_vegas_chisq(state));
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
  std::cout << "[Vegas::Generate] [DEBUG] " << fInputParameters->maxgen << " events will be generated" << std::endl;
  i = 0;
  while (i<fInputParameters->maxgen) {
    if (this->GenerateOneEvent()) i++;
  }
  std::cout << "[Vegas::Generate] [DEBUG] " << i << " events generated" << std::endl;
}

bool
Vegas::GenerateOneEvent()
{
  // Inherited from GMUGNA
  double ami, max;
  double y;
  int jj, jjj;
  double x[fFunction->dim];
  
  if (!fGenerationPrepared) {
    this->SetGen();
    fGenerationPrepared = true;
  }

  y = -1.;
  ami = 1./fMbin;
  max = pow(fMbin, fFunction->dim);

  // Correction cycles are started
  if (fJ!=0) {
  line4:
#ifdef DEBUG
    std::cout << "[Vegas::GenerateOneEvent] [DEBUG] Correction cycles are started."
	      << "\n\tj = " << fJ
	      << "\n\tcorrec = " << fCorrec
	      << "\n\tcorre2 = " << fCorrec2
	      << std::endl;
#endif
    if (fCorrec<1.) {
      if ((double)rand()/RAND_MAX>=fCorrec) {
        goto line7; //FIXME need to remove these goto
      }
      fCorrec = -1.;
    }
    else {
      fCorrec -= 1.;
    }
    // Select x values in Vegas bin
    for (unsigned int k=0; k<fFunction->dim; k++) {
      x[k] = ((double)rand()/RAND_MAX+fN[k])*ami;
    }
    // Compute weight for x value
    /*if (fInputParameters->ntreat>0) _weight = Treat(x);
    else _weight = this->F(x);*/
    _weight = F(x);
    // Parameter for correction of correction
    if (_weight>fFmax[fJ]) {
      if (_weight>_fmax2) _fmax2 = _weight;
      fCorrec2 -= 1.;
      fCorrec += 1.;
    }
    // Accept event
    if (_weight>=_fmdiff*(double)rand()/RAND_MAX+_fmold) { // FIXME!!!!
      return this->StoreEvent(x);
    }
    goto line4;
  line7:
    // Correction if too big weight is found while correction
    // (All your bases are belong to us...)
    if (_fmax2>fFmax[fJ]) {
      _fmold = fFmax[fJ];
      fFmax[fJ] = _fmax2;
      _fmdiff = _fmax2-_fmold;
      if (_fmax2<fFGlobalMax) {
        fCorrec = (_nm[fJ]-1.)*_fmdiff/fFGlobalMax-fCorrec2;
      }
      else {
        fFGlobalMax = _fmax2;
        fCorrec = (_nm[fJ]-1.)*_fmdiff/fFGlobalMax*_fmax2/fFGlobalMax-fCorrec2;
      }
      fCorrec2 = 0.;
      _fmax2 = 0.;
      goto line4;
      //return this->GenerateOneEvent(); //GOTO 4
    }
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
      x[i] = ((double)rand()/RAND_MAX+fN[i])*ami;
      jj = jjj;
    }
    
    // Get weight for selected x value
    /*if (fInputParameters->ntreat>0) _weight = this->Treat(x);
    else _weight = this->F(x);*/
    _weight = F(x);
    
    // Eject if weight is too low
    //if (y>_weight) {
    //std::cout << "ERROR : y>weight => " << y << ">" << _weight << ", " << fJ << std::endl;
      //_force_correction = false;
      //_force_returnto1 = true;
      //return this->GenerateOneEvent();
      //goto line1;
    //}
  } while (y>_weight);

  if (_weight<=fFmax[fJ]) fJ = 0;
  // Init correction cycle if weight is higher than fmax or ffmax
  else if (_weight<=fFGlobalMax) {
    _fmold = fFmax[fJ];
    fFmax[fJ] = _weight;
    _fmdiff = _weight-_fmold;
    fCorrec = (_nm[fJ]-1.)*_fmdiff/fFGlobalMax-1.;
  }
  else {
    _fmold = fFmax[fJ];
    fFmax[fJ] = _weight;
    _fmdiff = _weight-_fmold;
    fFGlobalMax = _weight;
    fCorrec = (_nm[fJ]-1.)*_fmdiff/fFGlobalMax*_weight/fFGlobalMax-1.;
  }
#ifdef DEBUG
  std::cout << "[Vegas::GenerateOneEvent] [DEBUG] correc = " << fCorrec << ", j = " << fJ << std::endl;
#endif
  // Return with an accepted event
  return this->StoreEvent(x);
}

bool
Vegas::StoreEvent(double *x_)
{
  if (_weight<=0.) {
#ifdef DEBUG
    std::cout << "[Vegas::StoreEvent] [DEBUG] Tried to store event while the weight is <= 0 : " << _weight << std::endl;
#endif
    return false;
  }
  fInputParameters->store = true;
  /*if (fInputParameters->ntreat>0) _weight = Treat(x_);
  else _weight = this->F(x_);*/
  _weight = F(x_);
  fInputParameters->ngen += 1;
  fInputParameters->store = false;
#ifdef DEBUG
  if (fInputParameters->ngen%1000==0) {
    std::cout << "[Vegas::StoreEvent] Generated events : " << fInputParameters->ngen << std::endl;
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
      /*if (fInputParameters->ntreat>0) z = this->Treat(x);
      else z = this->F(x);*/
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
    std::cout << "[Vegas::SetGen] [DEBUG] in iteration #" << i << " :"
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
  std::cout << "[Vegas::SetGen] [DEBUG]"
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

void
Vegas::DumpGrid()
{
  unsigned int i,j;
  // DUMP THE GRID
  for (i=0; i<fFunction->dim; i++) {
    for (j=0; j<fMaxNbins; j++) {
      std::cout << i << "\t" << j << "\t" << fCoord[j][i] << std::endl;
    }
  }
}

double
Vegas::Treat(double *x_, Parameters* ip_, bool storedbg_)
{
  double w, xx, y, dd, f;
  unsigned int i;
  int j;
  double z[fFunction->dim];

  if (_nTreatCalls==0) {
    _nTreatCalls = 1;
    _rTreat = std::pow(fStateBins, fFunction->dim);
    if (storedbg_ && remove("test_vegas")!=0) {
      std::cerr << "Error while trying to delete test_vegas" << std::endl;
    }
    //this->DumpGrid();
  }

  w = _rTreat;
  for (i=0; i<fFunction->dim; i++) {
    xx = x_[i]*fStateBins-1;
    j = xx;
    y = xx-j;
    if (j<=0) {
      dd = fCoord[0][i];
    }
    else {
      dd = fCoord[j+1][i]-fCoord[j][i];
    }
    z[i] = fCoord[j+1][i]-dd*(1.-y);
    w = w*dd;
  }

  f = this->F(z, ip_);

  if (storedbg_) {
    std::ofstream df;
    df.open("test_vegas", std::ios::app);
    df << w 
       << "\t" << w*f;
    for (unsigned int i=0; i<fFunction->dim; i++) {
      df << "\t" << z[i];
    }
    for (unsigned int i=0; i<fFunction->dim; i++) {
      df << "\t" << x_[i];
    }
    df << std::endl;
    df.close();
  }
#ifdef DEBUG
  std::cout << "[Vegas::Treat] [DEBUG] w = " << w << ", dd = " << dd << ", ndo = " << fStateBins << ", r = " << _rTreat << std::endl;
#endif
  return w*f;
}
