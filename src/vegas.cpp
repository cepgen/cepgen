#include "../include/vegas.h"

Vegas::Vegas(int dim_, double f_(double*,size_t,void*), InputParameters* inParam_)
{
  gsl_monte_vegas_params par;
  /* x content :
      0 = t1 mapping
      1 = t2 mapping
      2 = s2 mapping
      3 = yy4 definition
      4 = w4 mapping
      5 = xx6 definition
    ( 6 = phicm6 definition ) <- single-, double-dissociative only
    ( 7 = xq, wx mappings   ) <- double-dissociative only
  */

  _xl = new double[dim_];
  _xu = new double[dim_];
  
  _nTreatCalls = 0;

  for (int i=0; i<dim_; i++) {
    _xl[i] = 0.;
    _xu[i] = 1.;
  }
  _ndim = (size_t)dim_;
  _ncalls = 1e3;
  _nIter = 10;
  _s = gsl_monte_vegas_alloc(_ndim);
  gsl_monte_vegas_params_get(_s, &par);
  par.stage = 0;
  par.iterations = _nIter;
  par.verbose = -1;
  gsl_monte_vegas_params_set(_s, &par);
  gsl_rng_env_setup(); //FIXME ???
  _r = gsl_rng_alloc(gsl_rng_default);
  //gsl_rng_set(_r, 3001187); // sets the MC generator's seed

#ifdef DEBUG
  std::cout << "[Vegas::Vegas] [DEBUG]"
            << "\n  number of integration dimensions : " << dim_
            << "\n  number of iterations : " << _ncalls
            << std::endl;
#endif

  // GSL function to integrate on the whole phase space, provided with its input
  // parameters
  _F = new gsl_monte_function;
  _F->f=f_;
  _F->dim = (size_t)_ndim;
  _F->params = (void*)inParam_;

}

Vegas::~Vegas()
{
  gsl_monte_vegas_free(_s);
  gsl_rng_free(_r);

  delete[] _xl;
  delete[] _xu;
  delete _F;
}

int Vegas::Integrate(double *result_, double *abserr_)
{
  int vegas_status;
  gsl_monte_vegas_params par;

#ifdef DEBUG
  std::cout << "[Vegas::Integrate] [DEBUG] Vegas warm-up started !" << std::endl;
#endif

  // VEGAS warm up!
  gsl_monte_vegas_params_get(_s, &par);
  par.stage = 0;
  gsl_monte_vegas_params_set(_s, &par);
  vegas_status = gsl_monte_vegas_integrate(_F, _xl, _xu, _ndim, 1e3, _r, _s, result_, abserr_);

#ifdef DEBUG
  std::cout << "[Vegas::Integrate] [DEBUG] Vegas warm-up finished !" << std::endl;
#endif
  
  //FILE *fveg = fopen("vegas","w");
  
  gsl_monte_vegas_params_get(_s, &par);
  par.stage = 1;
  //par.iterations = _nIter;
  par.verbose = -1;
  //par.ostream = fveg;
  gsl_monte_vegas_params_set(_s, &par);

  //do {
    vegas_status = gsl_monte_vegas_integrate(_F, _xl, _xu, _ndim, _ncalls, _r, _s, result_, abserr_);
#ifdef DEBUG
    std::cout << "[Vegas::Integrate] [DEBUG] chisq/ndof = " << gsl_monte_vegas_chisq(_s) << std::endl;
#endif
  //} while (fabs(gsl_monte_vegas_chisq(_s)-1.)>.5);

#ifdef DEBUG
  std::cout << "[Vegas::Integrate] [DEBUG] (" << vegas_status << ") -> Computed cross-section = (" << *result_ << " +/- " << *abserr_ << ") pb" << std::endl;
#endif
  gsl_monte_vegas_params_get(_s, &par);
  //std::cout << "--> " << _s->bins << std::endl;
  return vegas_status;
}

void Vegas::SetGen()
{
  double sum, sum2, sum2p;
  double fsum, fsum2;
  double max, z;
  int n[_ndim];
  int i, j, jj, jjj, k, m;
  int npoin;
  double x[_ndim];
  double nm[7000], fmax[7000];
  double ffmax;
  double av, av2, sig, sig2, sigp, eff, eff1, eff2;
  
  npoin = 100;
  
  _mbin = 3.;
  ffmax = 0.;
  sum = sum2 = sum2p = 0.;
  max = pow(_mbin, _ndim);
  for (i=0; i<max; i++) {
    nm[i] = 0;
    fmax[i] = 0;
  }
  for (j=0; j<max; j++) {
    jj = j-1;
    for (k=0; k<(int)_ndim; k++) {
      jjj = jj/_mbin;
      n[k] = jj-jjj*_mbin;
      jj = jjj;
    }
    fsum = fsum2 = 0.;
    for (m=0; m<npoin; m++) {
      for (k=0; k<(int)_ndim; k++) {
        x[k] = ((double)rand()/(double)RAND_MAX+n[k])/_mbin;
      }
      if (_nTreatCalls>0) {
        z = Treat(x);
      }
      else {
        z = _F->f(x, _F->dim, _F->params);
      }
      if (z>fmax[j]) {
        fmax[j] = z;
      }
      fsum += z;
      fsum2 += pow(z, 2);
      /// XXXXXX.............
    }
    av = fsum/npoin;
    av2 = fsum2/npoin;
    sig2 = av2-pow(av, 2);
    sig = sqrt(sig2);
    sum += av;
    sum2 += av2;
    sum2p += sig2;
    if (fmax[j]>ffmax) {
      ffmax = fmax[j];
    }
    eff = 1.e4;
    if (fmax[j]!=0.) {
      eff = fmax[j]/av;
    }
  }
  sum /= max;
  sum2 /= max;
  sum2p /= max;
  sig = sqrt(sum2-pow(sum, 2));
  sigp = sqrt(sum2p);
  eff1 = 0.;
  for (i=0; i<max; i++) {
    eff1 += fmax[i];
  }
  eff1 /= (max*sum);
  eff2 = ffmax/sum;
  //std::cout << "ffmax = " << ffmax << std::endl;
  _ffmax = ffmax;
}

int Vegas::Generate(int nEvt_)
{
  int i;
  std::ofstream of;
  InputParameters *inp;
  GamGam *gg;
  
  SetGen();
  
  inp = (InputParameters*)_F->params;
  inp->store = false;
  inp->ngen = 0;
  inp->file = &of;
  _F->params = (void*)inp;
  of.open("test", std::ios_base::out);
  
  for (i=0; i<nEvt_; i++) {
    if (i%1000==0) {
      std::cout << "[Vegas::Generate] Event " << i << "/" << nEvt_ << std::endl;
    }
    if (!GenerateOneEvent(&of)) {
      continue;
    }
  }
  return 0;
}

//int Vegas::GenerateOneEvent(GamGam *gg)
int Vegas::GenerateOneEvent(std::ofstream *of)
{
  ((InputParameters*)_F->params)->store = false;
  // Inherited from GMUGNA
/*      INTEGER          NDIM,NCVG,ITMX,NPRN,IGRAPH,
     &                 NPOIN,NPRIN,NTREAT,IBEG,IEND
      COMMON /VEGPAR/  NDIM,NCVG,ITMX,NPRN,IGRAPH,
     &                 NPOIN,NPRIN,NTREAT,IBEG,IEND
      save /vegpar/*/

/**KEEP,COMGNA.
      INTEGER NGNA
      COMMON /COMMUP/  NGNA
      save /commup/*/
/*      COMMON/VGMAXI/MDUM,MBIN,FFMAX,FMAX(7000),NM(7000)
      save /vgmaxi/*/
  // VEGPAR
  int ntreat;
  // VGMAXI
  double fmax[7000], nm[7000];
  //
  double x[10], n[10];
  double weight, correc;
  double corre2, fmdiff, fmold, fmax2;
  int i, j, k;

  double ami, max;
  double y;
  double jj, jjj;
  //double ffmax;
  
  ntreat = 1;
  correc = 0.;
  
  //ffmax = Maximise();
  
  j = 0;
  ami = 1./_mbin;
  max = pow(_mbin, _ndim);
  
  // CORRECTION CYCLES ARE STARTED
  if (j!=0) {
line4:
    if (correc>=1.) {
      correc -= 1.;
    }
    if ((double)rand()/(double)RAND_MAX<correc) {
      correc = -1.;
      for (i=0; i<(int)_ndim; i++) {
        x[i] = ((double)rand()/(double)RAND_MAX+n[i])*ami;
      }
      if (ntreat>0) {
        weight = Treat(x);
      }
      else {
        weight = _F->f(x, _F->dim, _F->params);
      }
      if (weight>fmax[j]) {
        if (weight>fmax2) {
          fmax2 = weight;
        }
        corre2 -= 1.;
        correc += 1.;
      }
      if (weight>=fmdiff*(double)rand()/(double)RAND_MAX+fmold) {
//!             print *,'gmugna: done',weight,' x',(x(ipr),ipr=1,NDIM)
        return -1;
      }
      goto line4;
    }
    else { // 7
      if (fmax2>fmax[j]) {
        fmold = fmax[j];
        fmax[j] = fmax2;
        fmdiff = fmax2-fmold;
        if (fmax2<=_ffmax) {
          correc = (nm[j]-1.)*fmdiff/_ffmax-corre2;
        }
        else {
          _ffmax = fmax2;
          correc = (nm[j]-1.)*fmdiff/_ffmax*fmax2/_ffmax-corre2;
        }
        corre2 = 0.;
        fmax2 = 0.;
        goto line4;
      }
    }
  }
  // NORMAL GENERATION CYCLE STARTS HERE
  //std::cout << "max = " << max << std::endl;
  do {
    do {
      j = (int)(((double)rand()/(double)RAND_MAX)*max+1);
      y = (double)rand()/(double)RAND_MAX*_ffmax;
      nm[j] += 1;
    } while (y>fmax[j]);
    
    // SEL X VALUES IN THIS VEGAS BIN
    jj = j-1;
    for (k=0; k<(int)_ndim; k++) {
      jjj = jj/_mbin;
      n[k] = jj-jjj*_mbin;
      x[k] = ((double)rand()/(double)RAND_MAX+n[k])*ami;
      //std::cout << "x[" << k << "] = " << x[k] << std::endl;
      jj = jjj;
    }
    // GET WEIGHT FOR SEL X VALUES
    if (ntreat>0) {
      weight = Treat(x);
    }
    else {
      weight = _F->f(x, _F->dim, _F->params);
    }
    //InputParameters* ip = (InputParameters*)(_F->params);
    // REJECT IF WEIGHT IS TOO LOW
  } while (y>weight);
  if (weight<=fmax[j]) {
    j = 0;
  }
  // INIT CORRECTION CYCLES IF WEIGHT IS HIGHER THEN FMAX OR FFMAX
  else if (weight<=_ffmax) {
    fmold = fmax[j];
    fmax[j] = weight;
    fmdiff = weight-fmold;
    correc = (nm[j]-1.)*fmdiff/_ffmax-1.;
  }
  else {
    fmold = fmax[j];
    fmax[j] = weight;
    fmdiff = weight-fmold;
    _ffmax = weight;
    correc = (nm[j]-1.)*fmdiff/_ffmax*weight/_ffmax-1.;
  }
  InputParameters *inp = (InputParameters*)_F->params;
  inp->store = true;
  inp->file = of;
  //std::cout << "SetGen :: f called ! ngen = " << inp->ngen << std::endl;
  //std::cout << " --> f:params : store = " << ((InputParameters*)_F->params)->store << std::endl;
  //_F->params = (void*)inp;
  return _F->f(x, _F->dim, inp);
  //return 0;
}

double Vegas::Treat(double x[])
{
  int i;
  double w;
  double z[_ndim];
  int xx, jj;
  double y, dd;
  int j;
  if (_nTreatCalls==0) {
    _nTreatCalls = 1;
    _rTreat = pow(_s->bins, _ndim);
  }
  w = _rTreat;
  for (i=0; i<(int)_ndim; i++) {
    xx = x[i]*_s->bins;
    j = xx;
    jj = j+1;
    y = xx-j;
    if (j<=0.) {
      dd = _s->xi[0*_s->dim+i];
    }
    else {
      dd = _s->xi[jj*_s->dim+i]-_s->xi[j*_s->dim+i];
    }
    z[i] = _s->xi[jj*_s->dim+i]-dd*(1.-y);
    w *= dd;
  }
#ifdef DEBUG
  std::cout << "[Vegas::Treat] _s->bins = " << _s->bins << std::endl;
  std::cout << "[Vegas::Treat] dd at (i = " << i << " ; y = " << y << ") = " << dd << std::endl;
  std::cout << "[Vegas::Treat] f(z) = " << _F->f(z, _F->dim, _F->params) << std::endl;
  std::cout << "[Vegas::Treat] w = " << w << std::endl;
#endif
  return w*_F->f(z, _F->dim, _F->params);
}

int Vegas::LaunchGeneration(int nEvts_)
{
  int vegas_status;
  int niter, ngen;
  double result, abserr;
  std::ofstream of;
  gsl_monte_vegas_params par;
  
  std::cout << "[Vegas::LaunchGeneration] [DEBUG] Number of events to generate = " << nEvts_ << std::endl;

  InputParameters *inp = (InputParameters*)_F->params;
  inp->store = true;
  inp->ngen = 0;
  inp->file = &of;
  _F->params = (void*)inp;

  gsl_monte_vegas_params_get(_s, &par);
  par.stage = 3;
  par.iterations = 5;
  //par.mode = GSL_VEGAS_MODE_IMPORTANCE;
  par.verbose = -1;
  gsl_monte_vegas_params_set(_s, &par);

  of.open("test", std::ios_base::out);
#ifdef DEBUG
  gsl_monte_vegas_params_get(_s, &par);
  std::cout << "[Vegas::LaunchGeneration] [DEBUG] VEGAS stage after warm-up and cross-section computation" << par.stage << std::endl;
#endif

  niter = 0;
  ngen = 0;
  vegas_status = 0;
  while (ngen<nEvts_) {
    vegas_status = gsl_monte_vegas_integrate(_F, _xl, _xu, _ndim, (size_t)(nEvts_/(10*par.iterations)), _r, _s, &result, &abserr);
    ngen = ((InputParameters*)_F->params)->ngen;
#ifdef DEBUG
    std::cout << "[Vegas::LaunchGeneration] [DEBUG] Iteration number " << niter << std::endl;
    std::cout << "\t" << ngen << " events generated" << std::endl;
#endif
    niter++;
  }
  std::cout << "--> " << ngen << " events generated after " << niter << " iterations" << std::endl;
  
  of.close();
  /*std::cout << "[Vegas::LaunchGeneration] [DEBUG] chisq/ndof = " << gsl_monte_vegas_chisq(_s) << std::endl;
  std::cout << "[Vegas::LaunchGeneration] [DEBUG] (" << vegas_status << ") -> result : " << result << " +/- " << abserr << std::endl;*/
  
  return vegas_status;
}

