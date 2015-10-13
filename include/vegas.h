#ifndef _VEGAS_H
#define _VEGAS_H

#include <fstream>
#include <cstdio> // remove (DEBUG)
#include <gsl/gsl_monte_vegas.h>

#include "parameters.h"

#define fMaxNbins 50
#define ONE 1.

/**
 * Main occurence of the Monte-Carlo integrator@cite PeterLepage1978192 developed by G.P. Lepage in 1978
 * @brief Vegas Monte-Carlo integrator instance
 */
class Vegas {
  public:
    /**
     * Constructs the class by booking the memory and structures for the Vegas integrator. This code is based on the Vegas Monte Carlo integration algorithm developed by P. Lepage, as documented in @cite PeterLepage1978192
     * @param[in] dim_ The number of dimensions on which the function will be integrated
     * @param[in] f_ The function one is required to integrate
     * @param[inout] inParam_ A list of parameters to define the phase space on which this integration is performed (embedded in an Parameters object)
     */
    Vegas(const int dim_,double f_(double*,size_t,void*),Parameters* inParam_);
    /**
     * @brief Class destructor
     */
    ~Vegas();
    /**
     * Vegas algorithm to perform the n-dimensional Monte Carlo integration of a given function as described in @cite PeterLepage1978192
     * @author Primary author : G.P. Lepage
     * @author This C++ implementation : GSL
     * @param[out] result_ The cross section as integrated by Vegas for the given phase space restrictions
     * @param[out] abserr_ The error associated to the computed cross section
     * @return 0 if the integration was performed successfully
     */
    int Integrate(double* result_,double* abserr_);
    /**
     * Launches the Vegas generation of events according to the provided input
     *  parameters.
     * @brief Launches the generation of events
     */
    void Generate();
    /**
     * Generates one event according to the grid parameters set in Vegas::SetGen
     * @brief Generates one single event according to the method defined in the Fortran 77 version of LPAIR
     * @return A boolean stating if the generation was successful (in term of the computed weight for the phase space point)
     */
    bool GenerateOneEvent();
  private:
    /**
     * Evaluates the function to be integrated at a point @a x_, using the default Parameters object @a fInputParameters
     * @param[in] x_ The point at which the function is to be evaluated
     * @return Function value at this point @a x_
     */
    inline double F(double* x_) { return fFunction->f(x_, fFunction->dim, (void*)fInputParameters); }
    /**
     * Evaluates the function to be integrated at a point @a x_, given a set of Parameters @a ip_
     * @param[in] x_ The point at which the function is to be evaluated
     * @param[in] ip_ A set of parameters to fully define the function
     * @return Function value at this point @a x_
     */
    inline double F(double* x_,Parameters* ip_) { return fFunction->f(x_, fFunction->dim, (void*)ip_); }
    /**
     * Stores the event characterized by its _ndim-dimensional point in the phase
     * space to the output file
     * @brief Stores the event in the output file
     * @param[in] x_ The @a _ndim-dimensional point in the phase space defining the unique
     * event to store
     * @return A boolean stating whether or not the event could be saved
     */
    bool StoreEvent(double* x_);
    bool CorrectionCycle();
    /**
     * Sets all the generation mode variables and align them to the integration 
     * grid set while computing the cross-section
     * @brief Prepare the class for events generation
     */
    void SetGen();
    /**
     * @brief Integration grid size parameter
     */
    double fMbin;
    /**
     * @brief Lower bounds for the points to generate
     */
    double *fXlow;
    /**
     * @brief Upper bounds for the points to generate
     */
    double *fXup;
    /**
     * @brief Selected bin at which the function will be evaluated
     */
    int fJ;
    double fCorrec;
    double fCorrec2;
    /**
     * @brief List of parameters to specify the integration range and the physics determining the phase space
     */
    Parameters *fInputParameters;
    /**
     * @brief Flag to define whether or not the grid has been prepared for integration
     */
    bool fGridPrepared;
    /**
     * @brief Flag to define whether or not the generation has been prepared using @a SetGen (very time-consuming operation, thus needs to be called once)
     */
    bool fGenerationPrepared;
    /**
     * @brief Maximal value of the function at one given point
     */
    double *fFmax;
    /**
     * @brief Maximal value of the function in the considered integration range
     */
    double fFGlobalMax;
    int *fN;
    int *_nm;
    double _vegas_result;
    double _vegas_abserr;
    /**
     * @brief The function which will be integrated by this Vegas instance
     * @param x_ The point at which this function is evaluated
     * @param ndim_ The number of degrees of freedom this function has
     * @param params_ A "_void_-ified" Parameters object to define the boundaries of the phase space (physics constraints)
     */
    gsl_monte_function *fFunction;
    int fNumConverg;
    unsigned int fNumIter;
};

#endif

