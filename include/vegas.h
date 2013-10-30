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
    Vegas(int dim_,double f_(double*,size_t,void*),InputParameters* inParam_);
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
    int Integrate(double* result_,double* abserr_);
    /**
     * Launches the Vegas generation of events according to the provided input
     *  parameters.
     * @brief Launches the generation of events
     */
    int LaunchGeneration();
    void LaunchMyGeneration();
  private:
    //double Treat(double f_(double*,size_t,void*));
    //double Treat(gsl_monte_function*,double*);
    double Treat(double* x_,InputParameters* ip_);
    inline double Treat(double* x_) { return this->Treat(x_,(InputParameters*)(this->_F->params)); };
    //double Treat(double[]);
    /**
     * Stores the event characterized by its _ndim-dimensional point in the phase
     * space to the output file
     * @brief Stores the event in the output file
     * @brief x_ The _ndim-dimensional point in the phase space defining the unique
     * event to store
     * @return A boolean stating whether or not the event could be saved
     */
    bool StoreEvent(double*);
    /**
     * Generates one event according to the grid parameters set in Vegas::SetGen
     * @brief Generates one single event according to the method defined in the
     * Fortran 77 version of LPAIR
     * @return A boolean stating if the generation was successful (in term of the
     * computed weight for the phase space point)
     */
    bool GenerateOneEvent();
    //int GenerateOneEvent(GamGam*);
    //int GenerateOneEvent(std::ofstream*);
    /** @brief Number of times the Vegas::Treat method has been called */
    int _nTreatCalls;
    int _nTreat;
    double _rTreat;
    int _mbin;
    int *_n;
    int *_nm;
    double *_fmax;
    /**
     * Sets all the generation mode variables and align them to the integration 
     * grid set while computing the cross-section
     * @brief Prepare the class for events generation
     * @param of_ The file stream where to store the events after their generation
     */
    void SetGen(std::ofstream* of_);
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
    double _correc;
    double _weight;
    double _corre2;
    double _fmax2;
    double _fmdiff;
    double _fmold;
    int _j;
    InputParameters *_ip;
};

#endif

