#include "vegas.h"

#define COORD(s,i,j) ((s)->xi[(i)*(s)->dim + (j)])

Vegas::Vegas(int dim_, double f_(double*,size_t,void*), InputParameters* inParam_) :
  _nTreatCalls(0), _mbin(3),
  _ffmax(0.), _correc(0.), _corre2(0.), _fmax2(0.), _fmdiff(0.), _fmold(0.)
{
  gsl_monte_vegas_params par;
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

  _xl = new double[dim_];
  _xu = new double[dim_];

  // ...
  _n = new int[dim_];
  _nm = new int[7000];
  _fmax = new double[7000];
  // ...
  
  for (int i=0; i<dim_; i++) {
    _xl[i] = 0.;
    _xu[i] = 1.;
  }
  _ndim = (size_t)dim_;
  _nIter = inParam_->itvg;
  _ncalls = inParam_->ncvg;
  _s = gsl_monte_vegas_alloc(_ndim);
  gsl_monte_vegas_params_get(_s, &par);
  par.stage = 0;
  par.iterations = _nIter;
  par.verbose = -1;
  //par.mode = GSL_VEGAS_MODE_IMPORTANCE;
  //par.mode = GSL_VEGAS_MODE_STRATIFIED;
  //par.mode = GSL_VEGAS_MODE_IMPORTANCE_ONLY;
  //par.alpha = 1.;
  gsl_monte_vegas_params_set(_s, &par);
  gsl_rng_env_setup(); //FIXME ???
  _r = gsl_rng_alloc(gsl_rng_default);
  //gsl_rng_set(_r, 301187); // sets the MC generator's seed
  //
  // ...

#ifdef DEBUG
  std::cout << "[Vegas::Vegas] [DEBUG]"
            << "\n  Number of integration dimensions : " << dim_
            << "\n  Number of iterations : " << _nIter
            << "\n  Number of function calls : " << _ncalls
            << std::endl;
#endif

  // GSL function to integrate on the whole phase space, provided with its input
  // parameters
  _F = new gsl_monte_function;
  _F->f=f_;
  _F->dim = (size_t)_ndim;
  _F->params = (void*)inParam_;
  _ip = (InputParameters*)_F->params;
}

Vegas::~Vegas()
{
  gsl_monte_vegas_free(_s);
  gsl_rng_free(_r);

  delete _F;
  delete[] _xl;
  delete[] _xu;
  delete[] _n;
  delete[] _nm;
  delete[] _fmax;
}

int
Vegas::Integrate(double *result_, double *abserr_)
{
  int vegas_status;
  gsl_monte_vegas_params par;

#ifdef DEBUG
  std::cout << "[Vegas::Integrate] [DEBUG] Vegas warm-up started !" << std::endl;
#endif

  vegas_status = gsl_monte_vegas_init(_s);
  // VEGAS warm up!
  gsl_monte_vegas_params_get(_s, &par);
  par.stage = 0;
  par.verbose = -1;
  par.iterations = _ip->itvg;
  gsl_monte_vegas_params_set(_s, &par);
  //_ip->Dump();
  vegas_status = gsl_monte_vegas_integrate(_F, _xl, _xu, _ndim, 10000, _r, _s, result_, abserr_);

#ifdef DEBUG
  std::cout << "[Vegas::Integrate] [DEBUG] Vegas warm-up finished !" << std::endl;
#endif
  
  //FILE *fveg = fopen("vegas","w");
  gsl_monte_vegas_params_get(_s, &par);
  par.stage = 1;
  par.verbose = 0;
  //par.ostream = fveg;
  gsl_monte_vegas_params_set(_s, &par);

  //do {
  vegas_status = gsl_monte_vegas_integrate(_F, _xl, _xu, _ndim, _ncalls, _r, _s, result_, abserr_);
#ifdef DEBUG
  std::cout << "[Vegas::Integrate] [DEBUG] chisq/ndof = " << gsl_monte_vegas_chisq(_s) << std::endl;
#endif
  //} while (fabs(gsl_monte_vegas_chisq(_s)-1.)>.1);

#ifdef DEBUG
  std::cout << "[Vegas::Integrate] [DEBUG] (" << vegas_status << ") -> Computed cross-section = (" 
            << *result_ << " +/- "
            << *abserr_ << ") pb"
            << std::endl;
#endif
  /*gsl_monte_vegas_params_get(_s, &par);
  std::cout << "--> " << _s->bins << std::endl;
  std::cout << "--> " << _s->boxes << std::endl;*/
  return vegas_status;
}

int
Vegas::LaunchGeneration()
{
  int vegas_status;
  int niter, ngen;
  double result, abserr;
  std::ofstream of, fd;
  gsl_monte_vegas_params par; 

  std::cout << "[Vegas::LaunchGeneration] [DEBUG] Number of events to generate = " << _ip->maxgen << std::endl;

  _ip->store = true;
  _ip->ngen = 0;
  _ip->file = &of;
  //_ip->file_debug = &fd;

  gsl_monte_vegas_params_get(_s, &par);
  par.stage = 2;
  par.iterations = 1;
  //par.mode = GSL_VEGAS_MODE_IMPORTANCE;
  par.verbose = -1;
  gsl_monte_vegas_params_set(_s, &par);

  of.open("test", std::ios_base::out);
  //fd.open("debug", std::ios_base::out);
#ifdef DEBUG
  gsl_monte_vegas_params_get(_s, &par);
  std::cout << "[Vegas::LaunchGeneration] [DEBUG] VEGAS stage after warm-up and cross-section computation : " << par.stage << std::endl;
#endif
  std::cout << "dump after xsec computation" << std::endl;
  _ip->Dump();

  niter = 0;
  ngen = 0;
  vegas_status = 0;
  while (ngen<_ip->maxgen) {
    //vegas_status = gsl_monte_vegas_integrate(_F, _xl, _xu, _ndim, (size_t)(_ip->maxgen/(10*par.iterations)), _r, _s, &result, &abserr);
    vegas_status = gsl_monte_vegas_integrate(_F, _xl, _xu, _ndim, (size_t)_ip->maxgen, _r, _s, &result, &abserr);
    ngen = _ip->ngen;
//#ifdef DEBUG
    std::cout << "[Vegas::LaunchGeneration] [DEBUG] Iteration number " << niter
              << "\n\t" << ngen << " events generated"
              << std::endl;
//#endif
    niter++;
  }
  std::cout << "--> " << ngen << " events generated after " << niter << " iteration(s)" << std::endl;
  
  of.close();
#ifdef DEBUG
  std::cout << "[Vegas::LaunchGeneration] [DEBUG] chisq/ndof = " << gsl_monte_vegas_chisq(_s) << std::endl;
  std::cout << "[Vegas::LaunchGeneration] [DEBUG] (" << vegas_status << ") -> result : " << result << " +/- " << abserr << std::endl;
#endif

  return vegas_status;
}

void
Vegas::LaunchMyGeneration()
{
  //count_ = 1;
  std::ofstream of;
  std::string fn;
  int i;
  
  fn = "test";
  of.open(fn.c_str());
  this->SetGen(&of);
  std::cout << "[Vegas::LaunchMyGeneration] [DEBUG] " << _ip->maxgen << " events will be generated" << std::endl;
  i = 0;
  while (i<_ip->maxgen) {
    if (this->GenerateOneEvent()) i++;
  }
  std::cout << "[Vegas::LaunchMyGeneration] [DEBUG] " << i << " events generated in \"" << fn << "\"" << std::endl;
  of.close();
}

bool
Vegas::GenerateOneEvent()
{
  // ...
  // Inherited from GMUGNA
  double ami, max;
  double y;
  int jj, jjj;
  double x[_ndim];
  
  ami = 1./_mbin;
  max = pow(_mbin, _ndim);
  //std::cout << "ffmax = " << _ffmax << ", fmax2 = " << _fmax2 << ", j = " << _j << ", fmold = " << _fmold << std::endl;

  // Correction cycles are started
  if (_j!=0) {
#ifdef DEBUG
    std::cout << "[Vegas::GenerateOneEvent] [DEBUG] Correction cycles are started."
	      << "\n\tj = " << _j
	      << "\n\tcorrec = " << _correc
	      << "\n\tcorre2 = " << _corre2
	      << std::endl;
#endif
    if (_correc<1.) {
      if ((double)rand()/RAND_MAX<_correc) { //GOTO 7
	std::cout << ">>>> GOTO 7" << std::endl;
        // Correction if too big weight is found while correction
        // (All your bases are belong to us...)
        if (_fmax2>_fmax[_j]) {
          _fmold = _fmax[_j];
          _fmax[_j] = _fmax2;
          _fmdiff = _fmax2-_fmold;
          if (_fmax2<_ffmax) {
            _correc = (_nm[_j]-1.)*_fmdiff/_ffmax-_corre2;
          }
          else {
            _ffmax = _fmax2;
            _correc = (_nm[_j]-1.)*_fmdiff/_ffmax*_fmax2/_ffmax-_corre2;
          }
	  std::cout << "new correc = " << _correc << std::endl;
          _corre2 = 0.;
          _fmax2 = 0.;
          //std::cout << "_j => " << _j << std::endl;
          //_j = -1; //FIXME to ensure the first condition is respected
          return this->GenerateOneEvent(); //GOTO 4
        }
      }
      _correc = -1.;
    }
    else {
      _correc -= 1.;
    }
    // Select x values in Vegas bin
    for (unsigned int k=0; k<_ndim; k++) {
      x[k] = ((double)rand()/RAND_MAX+_n[k])*ami;
    }
    // Compute weight for x value
    if (_ip->ntreat>0) _weight = Treat(x);
    else _weight = this->F(x);
    // Parameter for correction of correction
    if (_weight>_fmax[_j]) {
      if (_weight>_fmax2) _fmax2 = _weight;
      _corre2 -= 1.;
      _correc += 1.;
    }
    // Accept event
    if (_weight>=_fmdiff*(double)rand()/RAND_MAX+_fmold) { // FIXME!!!!
      //std::cout << "-------> yes !" << std::endl;
      return this->StoreEvent(x);
      //return true;
    }
    //std::cout << "----> no ! new event to be generated" << std::endl;
    return this->GenerateOneEvent();
  }
  // Normal generation cycle
  // Select a Vegas bin and reject if fmax is too little
  y = 1.e10;
  do {
    // ...
    _j = (double)rand()/RAND_MAX*max;
    y = (double)rand()/RAND_MAX*_ffmax;
    _nm[_j] += 1;
  } while (y>_fmax[_j]);
  //std::cout << y << std::endl;
  // Select x values in this Vegas bin
  //jj = _j-1;
  jj = _j; //FIXME!!!!!!! Ugly ugly bad bad
  for (unsigned int i=0; i<_ndim; i++) {
    jjj = jj/_mbin;
    _n[i] = jj-jjj*_mbin;
    x[i] = ((double)rand()/RAND_MAX+_n[i])*ami;
    //std::cout << "x[" << i << "] = " << x[i] << ", jj = " << jj << ", jjj = " << jjj << std::endl;
    jj = jjj;
  }
  //std::cout << "jjj = " << jjj << ", jj = " << jj << ", j-1 = " << _j-1 << std::endl;

  // Get weight for selected x value
  if (_ip->ntreat>0) _weight = Treat(x);
  else _weight = this->F(x);
  //std::cout << ">>> " << _weight << std::endl;

  // Eject if weight is too low
  if (y>_weight) {
    //std::cout << "ERROR : y>weight => " << y << ">" << _weight << std::endl;
    _j = 0;
    return this->GenerateOneEvent();
    //return false;
  }

  if (_weight<=_fmax[_j]) _j = 0;
  // Init correction cycle if weight is higher than fmax or ffmax
  else if (_weight<=_ffmax) {
    _fmold = _fmax[_j];
    _fmax[_j] = _weight;
    _fmdiff = _weight-_fmold;
    _correc = (_nm[_j]-1.)*_fmdiff/_ffmax-1.;
  }
  else {
    _fmold = _fmax[_j];
    _fmax[_j] = _weight;
    _fmdiff = _weight-_fmold;
    _ffmax = _weight;
    _correc = (_nm[_j]-1.)*_fmdiff/_ffmax*_weight/_ffmax-1.;
  }
  //#ifdef DEBUG
  std::cout << "[Vegas::GenerateOneEvent] [DEBUG] correc = " << _correc << ", j = " << _j << std::endl;
  //#endif
  std::cout << "accepted event !" << std::endl;
  // Return with an accepted event
  return this->StoreEvent(x);
  //return true;
}

bool
Vegas::StoreEvent(double *x_)
{
  if (_weight<=0.) {
//#ifdef DEBUG
    std::cout << "[Vegas::StoreEvent] [DEBUG] Tried to store event while the weight is <= 0 : " << _weight << std::endl;
//#endif
    return false;
  }
  _ip->store = true;
  //std::cout << "-------> storing the event !" << std::endl;
  if (_ip->ntreat>0) _weight = Treat(x_);
  else _weight = this->F(x_);
  //std::cout << _weight-_F->f(x_, _ndim, (void*)_F->params) << std::endl;
  //_F->f(x_, _ndim, (void*)_F->params);
  _ip->ngen += 1;
  _ip->store = false;
  if (_ip->ngen%10000==0) {
    std::cout << "[Vegas::StoreEvent] Generated events : " << _ip->ngen << std::endl;
  }
  //std::cout << "--> " << _ip->ngen << std::endl;

  return true;
}

void
Vegas::SetGen(std::ofstream *of_)
{
  int max;
  //int jj, jjj;
  int jj;
  double jjj;
  double sum, sum2, sum2p;
  int n[10];
  int npoin = _ip->npoints;
  double fsum, fsum2;
  double z;
  double x[_ndim];
  double sig2;
  double av, av2;
  //#ifdef DEBUG
  double eff, eff1, eff2;
  double sig, sigp;
  //#endif

  _ip->ngen = 0;
  _ip->file = of_;
#ifdef DEBUG
  std::cout << "[Vegas::SetGen] [DEBUG] maxgen = " << _ip->maxgen << std::endl;
  _ip->Dump();
#endif
  // ...
  sum = 0.;
  sum2 = 0.;
  sum2p = 0.;
  max = pow(_mbin, _ndim);

  for (int i=0; i<max; i++) {
    _nm[i] = 0;
    _fmax[i] = 0.;
  }

  for (int i=0; i<max; i++) {
    //jj = i-1;
    jj = i; //FIXME !!!!!
    for (unsigned int j=0; j<_ndim; j++) {
      //jjj = floor(jj/_mbin); //FIXME floor or ceil???
      jjj = jj/_mbin;
      n[j] = jj-jjj*_mbin;
      //std::cout << "--> j(" << j << ") : n = " << n[j] << ", jj = " << jj << ", jjj = " << jjj << std::endl;
      jj = jjj+1;
    }
    fsum = 0.;
    fsum2 = 0.;
    for (int j=0; j<npoin; j++) {
      for (unsigned int k=0; k<_ndim; k++) {
        x[k] = ((double)rand()/RAND_MAX+n[k])/_mbin;
      }
      if (_ip->ntreat>0) z = this->Treat(x);
      else z = this->F(x);
      if (z<0) std::cout << "z=" << z << std::endl;
      if (z>_fmax[i]) _fmax[i] = z;
      fsum += z;
      fsum2 += pow(z, 2);
      //std::cout << fsum << std::endl;
    }
    //std::cout << "fsum=" << fsum << std::endl;
    av = fsum/npoin;
    av2 = fsum2/npoin;
    sig2 = av2-pow(av, 2);
    //#ifdef DEBUG
    sig = sqrt(sig2);
    //#endif
    sum += av;
    sum2 += av2;
    sum2p += sig2;
    if (_fmax[i]>_ffmax) _ffmax = _fmax[i];
    //#ifdef DEBUG
    eff = 1.e4;
    if (_fmax[i]!=0.) eff = _fmax[i]/av;
    //std::cout << "fmax(" << i << ") = " << _fmax[i] << std::endl;
    //std::cout << "fsum = " << fsum << ", fsum2 = " << fsum2 << std::endl;
    /*std::cout << "[Vegas::SetGen] [DEBUG] in iteration #" << i << " :"
	      << "\n\tav   = " << av
	      << "\n\tsig  = " << sig
	      << "\n\tfmax = " << _fmax[i]
	      << "\n\teff  = " << eff
	      << "\n\tn = (";
    for (unsigned int j=0; j<_ndim; j++) {
      std::cout << n[j];
      if (j!=_ndim-1) std::cout << ", ";
    }
    std::cout << ")" << std::endl;*/
    //#endif
  }

  sum = sum/max;
  sum2 = sum2/max;
  sum2p = sum2p/max;
  //std::cout << "sum = " << sum << ", sum2 = " << sum2 << ", sum2p = " << sum2p << std::endl;

  //#ifdef DEBUG
  sig = sqrt(sum2-pow(sum, 2));
  sigp = sqrt(sum2p);
  eff1 = 0.;
  for (int i=0; i<max; i++) {
    eff1 += _fmax[i];
  }
  eff1 = eff1/(max*sum);
  eff2 = _ffmax/sum;
  std::cout << "[Vegas::SetGen] [DEBUG]"
            << "\n\tAverage function value     = sum   = " << sum
            << "\n\tOverall standard deviation = sig   = " << sig
	    << "\n\tAverage standard deviation = sigp  = " << sigp
	    << "\n\tMaximum function value     = ffmax = " << _ffmax
	    << "\n\tAverage inefficiency       = eff1  = " << eff1 
	    << "\n\tOverall inefficiency       = eff2  = " << eff2 
            << "\n\teff = " << eff 
	    << std::endl;
  //#endif
}

double
Vegas::Treat(double *x_, InputParameters* ip_)
{
  double w;
  double xx;
  int j, jj;
  int ndo = _s->bins;
  double y;
  double dd;
  double z[_ndim];

  if (_nTreatCalls==0) {
    _nTreatCalls = 1;
    _rTreat = pow(ndo, _ndim);
  }
  w = _rTreat;
  //std::cout << "rtreat = " << _rTreat << std::endl;
  for (unsigned int i=0; i<_ndim; i++) {
    xx = x_[i]*ndo;
    j = xx;
    jj = j+1;
    y = xx-j;
    //std::cout << "y[" << i << "] = " << y << ", jj = " << jj << ", ndo = " << ndo << std::endl;
    if (j<=0) {
      dd = COORD(_s,0,i);
    }
    else {
      //std::cout << "J>0" << std::endl;
      dd = COORD(_s,jj,i)-COORD(_s,j,i);
    }
    z[i] = COORD(_s,jj,i)-dd*(1.-y);
    w = w*dd;
    //std::cout << "dd[" << i << "]-> z[" << i << "]=" << z[i] << ", dd=" << dd << " --> w=" << w << ", y=" << y << ", j=" << j << ", jj=" << jj << std::endl;
  }
  //std::cout << "weight = " << w << std::endl;
  //ip_->Dump();
#ifdef DEBUG
  std::cout << "[Vegas::Treat] [DEBUG] w = " << w << ", dd = " << dd << ", ndo = " << ndo << ", r = " << _rTreat << std::endl;
#endif
  return w*this->F(z, ip_);
}

