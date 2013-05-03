#include "../include/vegas.h"

Vegas::Vegas(int dim_, double f_(double*,size_t,void*), InputParameters* inParam_)
{
  gsl_monte_vegas_params par;
  // x content : 
  //  0 : t1 mapping
  //  1 : t2 mapping
  //  2 : s2 mapping
  //  3 : yy4 definition
  //  4 : w4 mapping
  //  5 : xx6 definition
  // (6 : phicm6 definition)
  // (7 : xq, wx mappings)

  _xl = new double[dim_];
  _xu = new double[dim_];
  
  for (int i=0; i<dim_; i++) {
    _xl[i] = 0.;
    _xu[i] = 1.;
  }
  _ndim = (size_t)dim_;
  _ncalls = 5e3;
  _nIter = 10;
  _s = gsl_monte_vegas_alloc(_ndim);
  gsl_monte_vegas_params_get(_s, &par);
  par.stage = 0;
  par.iterations = _nIter;
  par.verbose = -1;
  gsl_monte_vegas_params_set(_s, &par);
  _r = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(_r, 3001187); // sets the MC generator's seed
  
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

int Vegas::LaunchIntegration()
{
  int vegas_status;
  double result, abserr;  
  gsl_monte_vegas_params par;
  
  std::cout << "Vegas warm-up started !" << std::endl;
  
  // VEGAS warm up!
  vegas_status = gsl_monte_vegas_integrate(_F, _xl, _xu, _ndim, 1e3, _r, _s, &result, &abserr);
  
  std::cout << "Vegas warm-up finished !" << std::endl;
  
  gsl_monte_vegas_params_get(_s, &par);
  par.stage = 1;
  //par.iterations = _nIter;
  par.verbose = 0;
  gsl_monte_vegas_params_set(_s, &par);
  
  //do {
    vegas_status = gsl_monte_vegas_integrate(_F, _xl, _xu, _ndim, _ncalls, _r, _s, &result, &abserr);
    std::cout << "[Vegas::LaunchIntegration] [DEBUG] chisq/ndof = " << gsl_monte_vegas_chisq(_s) << std::endl;
  //} while (fabs(gsl_monte_vegas_chisq(_s)-1.)>.5);
  
  std::cout << "[Vegas::LaunchIntegration] [DEBUG] (" << vegas_status << ") -> result : " << result << " +/- " << abserr << std::endl;
  
  return vegas_status;
}
int Vegas::LaunchGeneration(int nEvts_)
{
  int vegas_status;
  double result, abserr;
  std::ofstream of;
  gsl_monte_vegas_params par;

  std::cout << "[Vegas::LaunchGeneration] [DEBUG] Number of events to generate = " << nEvts_ << std::endl;

  InputParameters *inp = (InputParameters*)_F->params;
  inp->generation = true;
  inp->file = &of;
  _F->params = (void*)inp;
  
  gsl_monte_vegas_params_get(_s, &par);
  par.stage = 3;
  gsl_monte_vegas_params_set(_s, &par);
  
  of.open("test", std::ios_base::out);
#ifdef DEBUG
  gsl_monte_vegas_params_get(_s, &par);
  std::cout << "[Vegas::LaunchGeneration] [DEBUG] VEGAS stage after warm-up and cross-section computation" << par.stage << std::endl;
#endif
  
  vegas_status = gsl_monte_vegas_integrate(_F, _xl, _xu, _ndim, (size_t)(nEvts_/_nIter), _r, _s, &result, &abserr);
  /*std::cout << "[Vegas::LaunchGeneration] [DEBUG] chisq/ndof = " << gsl_monte_vegas_chisq(_s) << std::endl;
  std::cout << "[Vegas::LaunchGeneration] [DEBUG] (" << vegas_status << ") -> result : " << result << " +/- " << abserr << std::endl;*/
  return vegas_status;
}

