#ifndef _VEGAS_H
#define _VEGAS_H

#include <fstream>
#include <gsl/gsl_monte_vegas.h>

#include "gamgam.h"

/** @brief Vegas Monte-Carlo integrator instance */
class Vegas {
  public:
    /**
     * Constructs the class by booking the memory and structures for the
     *  GSL Vegas integrator. This code from the GNU scientific library is based
     *  on the Vegas Monte Carlo integration algorithm developed by P. Lepage.
     * @cite PeterLepage1978192
     * @param dim_ The number of dimensions on which the function will
     *  be integrated
     * @param f_ The function one is required to integrate
     * @param inParam_ A list of parameters to define the phase space on which
     *  this integration is performed (embedded in an InputParameters object)
     */
    Vegas(int,double f_(double*,size_t,void*),InputParameters* inParam_);
    /**
     * @brief Class destructor
     */
    ~Vegas();
    /**
     * Launches the Vegas integration of the provided function with the
     *  provided input parameters.
     * @brief Launches the integration of the provided function
     * @param result_ The cross section as integrated by Vegas for the given
     *  phase space restrictions
     * @param abserr_ The error associated to the computed cross section
     */
    int Integrate(double*,double*);
    /**
     * Launches the Vegas generation of events according to the provided input
     *  parameters.
     * @brief Launches the generation of events
     * @param nEvts_ The number of events to generate
     */
    int LaunchGeneration(int);
    int Generate(int);
    void SetGen();
  private:
    double Treat(double[]);
    //int GenerateOneEvent(GamGam*);
    int GenerateOneEvent(std::ofstream*);
    int _nTreatCalls;
    double _rTreat;
    double _mbin;
    double _ffmax;
    /** @brief GSL's random number generator */
    gsl_rng *_r;
    /** @brief GSL's Vegas integration state structure */
    gsl_monte_vegas_state *_s;
    /** @brief The wrapped-up function to integrate, along with the input parameters */
    gsl_monte_function *_F;
    /** @brief The number of dimensions on which to integrate the function */
    size_t _ndim;
    /** @brief Number of points to generate in order to integrate the function */
    size_t _nIter;
    /** @brief Fixed number of function calls to use */
    size_t _ncalls;
    /** @brief Lower bounds for the points to generate */
    double *_xl;
    /** @brief Upper bounds for the points to generate */
    double *_xu;
};

#endif

