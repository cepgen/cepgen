#include "vegas.h"

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
  _ncalls = 14000; // equivalent in LPAIR : NCVG->NCALLS
  _nIter = inParam_->itmx;
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
  //delete inp;

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

