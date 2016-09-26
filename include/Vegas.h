#ifndef Vegas_h
#define Vegas_h

#include <fstream>
#include <cstdio> // remove (DEBUG)
#include <gsl/gsl_monte_vegas.h>

#include "Parameters.h"

#define fMaxNbins 50
#define ONE 1.

/**
 * Main occurence of the Monte-Carlo integrator @cite PeterLepage1978192 developed by G.P. Lepage in 1978
 * \brief Vegas Monte-Carlo integrator instance
 */
class Vegas {
  public:
    /**
     * Book the memory slots and structures for the Vegas integrator
     * \note This code is based on the Vegas Monte Carlo integration algorithm developed by P. Lepage, as documented in @cite PeterLepage1978192
     * \param[in] dim_ Number of dimensions on which the function will be integrated
     * \param[in] f_ Function to be integrated
     * \param[inout] inParam_ Run parameters to define the phase space on which this integration is performed (embedded in an Parameters object)
     */
    Vegas(const int dim_,double f_(double*,size_t,void*),Parameters* inParam_);
    /// Class destructor
    ~Vegas();
    /**
     * Vegas algorithm to perform the n-dimensional Monte Carlo integration of a given function as described in @cite PeterLepage1978192
     * \author Primary author: G.P. Lepage
     * \author This C++ implementation: GSL
     * \param[out] result_ The cross section as integrated by Vegas for the given phase space restrictions
     * \param[out] abserr_ The error associated to the computed cross section
     * \return 0 if the integration was performed successfully
     */
    int Integrate(double* result_,double* abserr_);
    /// Launch the generation of events
    void Generate();
    /**
     * Generate one event according to the grid parameters set in Vegas::SetGen
     * \brief Generate one single event according to the method defined in the Fortran 77 version of LPAIR
     * \return A boolean stating if the generation was successful (in term of the computed weight for the phase space point)
     */
    bool GenerateOneEvent();
  private:
    /**
     * Evaluate the function to be integrated at a point @a x_, using the default Parameters object @a fInputParameters
     * \param[in] x_ The point at which the function is to be evaluated
     * \return Function value at this point @a x_
     */
    inline double F(double* x_) { return fFunction->f(x_, fFunction->dim, (void*)fInputParameters); }
    /**
     * Evaluate the function to be integrated at a point @a x_, given a set of Parameters @a ip_
     * \param[in] x_ The point at which the function is to be evaluated
     * \param[in] ip_ A set of parameters to fully define the function
     * \return Function value at this point @a x_
     */
    inline double F(double* x_,Parameters* ip_) { return fFunction->f(x_, fFunction->dim, (void*)ip_); }
    /**
     * Store the event characterized by its _ndim-dimensional point in the phase
     * space to the output file
     * @brief Store the event in the output file
     * @param[in] x_ The @a _ndim-dimensional point in the phase space defining the unique
     * event to store
     * @return A boolean stating whether or not the event could be saved
     */
    bool StoreEvent(double* x_);
    /// Start the correction cycle on the grid
    /// \param x_ Point in the phase space considered
    bool CorrectionCycle(double* x_);
    /**
     * Set all the generation mode variables and align them to the integration grid set while computing the cross-section
     * @brief Prepare the class for events generation
     */
    void SetGen();
    /// Integration grid size parameter
    double fMbin;
    /// Lower bounds for the points to generate
    double *fXlow;
    /// Upper bounds for the points to generate
    double *fXup;
    /// Selected bin at which the function will be evaluated
    int fJ;
    double fCorrec;
    double fCorrec2;
    /// List of parameters to specify the integration range and the physics determining the phase space
    Parameters *fInputParameters;
    /// Has the grid been prepared for integration?
    bool fGridPrepared;
    /// Has the generation been prepared using @a SetGen call? (very time-consuming operation, thus needs to be called once)
    bool fGenerationPrepared;
    bool fHasCorrection;
    /// Maximal value of the function at one given point
    double *fFmax;
    double fFmax2;
    double fFmaxDiff;
    double fFmaxOld;
    /// Maximal value of the function in the considered integration range
    double fFGlobalMax;
    int *fN;
    int *fNm;
    /// GSL structure storing the function to be integrated by this Vegas instance (along with its parameters)
    gsl_monte_function *fFunction;
    /// Number of function calls to be computed for each point
    int fNumConverg;
    /// Number of iterations for the integration
    unsigned int fNumIter;
    /// Generic array of components to be used in all parts of the integration and generation code
    double* fX;
};

#endif

