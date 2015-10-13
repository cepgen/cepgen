#ifndef _VEGAS_H
#define _VEGAS_H

#include <fstream>
#include <cstdio> // remove (DEBUG)

#include "parameters.h"

#define MAX_ND 50
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
     * Vegas algorithm to perform the (_dim)-dimensional Monte Carlo integration of a given function as described in @cite PeterLepage1978192
     * @author Primary author : G.P. Lepage
     * @author This C++ implementation : L. Forthomme
     * @date Sep 1976
     * @date Reviewed in Apr 1978
     * @date FTN5 version 21 Aug 1984
     * @date This C++ implementation is from 12 Dec 2013
     * @param[out] result_ The cross section as integrated by Vegas for the given phase space restrictions
     * @param[out] abserr_ The error associated to the computed cross section
     * @return 0 if the integration was performed successfully
     */
    int Integrate(double* result_,double* abserr_);
    /**
     * First stage of the integration process : Initialization of cumulative variables (no grid so far) 
     * @param ncalls_ Number of function calls to be performed
     * @return 0, if this part, along with parts 2 and 3 were performed successfully 
     */
    int Vegas1(int ncalls_=-1);
    /**
     * Second stage of the integration process : Grid initialization
     * @param ncalls_ Number of function calls to be performed
     * @return 0, if this part, along with part 3 were performed successfully 
     */
    int Vegas2(int ncalls_=-1);
    /**
     * Third stage of the integration process : Main loop 
     * @return 0, if this part was performed successfully 
     */
    int Vegas3();
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
     * Transforms the function to integrate into a numerically stable function where poles are tamed.
     * @param[in] x_ The @a _ndim -dimensional point at which the stabilised function is to be evaluated
     * @param[in] ip_ The physics parameters to apply to the function to evaluate
     * @param[in] storedbg_ A debugging flag to set whether or not the internal variables of this method need to be stored for further processing
     * @return Tamed function value at this point @a x_
     */
    double Treat(double* x_,Parameters* ip_,bool storedbg_=false);
    /**
     * Evaluates the smoothed version of the function to be integrated at a point @a x_, using the default Parameters object @a _ip 
     * @param[in] x_ The point at which the tamed function is to be evaluated
     * @return Tamed function value at this point @a x_
     */
    inline double Treat(double* x_) { return this->Treat(x_,(Parameters*)this->_ip); }
    /**
     * Evaluates the smoothed version of the function to be integrated at a point @a x_
     * @param[in] x_ The point at which the tamed function is to be evaluated
     * @param[in] storedbg_ A debugging flag to set whether or not the internal variables of this method need to be stored for further processing
     * @return Tamed function value at this point @a x_
     */
    inline double Treat(double* x_,bool storedbg_) { return this->Treat(x_,(Parameters*)this->_ip,storedbg_); }
    /**
     * Evaluates the function to be integrated at a point @a x_, using the default Parameters object @a _ip
     * @param[in] x_ The point at which the function is to be evaluated
     * @return Function value at this point @a x_
     */
    inline double F(double* x_) { return this->_f(x_, this->_ndim, (void*)this->_ip); }
    /**
     * Evaluates the function to be integrated at a point @a x_, given a set of Parameters @a ip_
     * @param[in] x_ The point at which the function is to be evaluated
     * @param[in] ip_ A set of parameters to fully define the function
     * @return Function value at this point @a x_
     */
    inline double F(double* x_,Parameters* ip_) { return this->_f(x_, this->_ndim, (void*)ip_); }
    /**
     * Stores the event characterized by its _ndim-dimensional point in the phase
     * space to the output file
     * @brief Stores the event in the output file
     * @param[in] x_ The @a _ndim-dimensional point in the phase space defining the unique
     * event to store
     * @return A boolean stating whether or not the event could be saved
     */
    bool StoreEvent(double* x_);
    /**
     * Sets all the generation mode variables and align them to the integration 
     * grid set while computing the cross-section
     * @brief Prepare the class for events generation
     */
    void SetGen();
    /**
     * Debugging method used to dump the integration grid in the standard output stream.
     */
    void DumpGrid();
    /**
     * @brief The number of dimensions on which to integrate the function
     */
    const size_t _ndim;
    unsigned int _ndo;
    /**
     * Has the Treat function already been called once ?
     */
    int _nTreatCalls;
    /**
     * \f$r = \text{ndo}^\text{ndim}\f$ value of the Treat function
     */
    double _rTreat;
    /**
     * @brief Integration grid size parameter
     */
    double _mbin;
    /**
     * @brief Maximal value of the function in the considered integration range
     */
    double _ffmax;
    int *_n;
    int *_nm;
    /**
     * @brief Maximal value of the function at one given point
     */
    double *_fmax;
    /**
     * @brief Lower bounds for the points to generate
     */
    double *_xl;
    /**
     * @brief Upper bounds for the points to generate
     */
    double *_xu;
    double _correc;
    /**
     * @brief Weight of the point in the total integration
     */
    double _weight;
    double _corre2;
    /**
     * @brief Square of the maximal function value in the integration grid
     * */
    double _fmax2;
    double _fmdiff;
    double _fmold;
    /**
     * @brief Selected bin at which the function will be evaluated
     */
    int _j;
    double *_xi[MAX_ND];
    double *_d[MAX_ND];
    double *_di[MAX_ND];
    /**
     * @brief List of parameters to specify the integration range and the physics determining the phase space
     */
    Parameters *_ip;
    /**
     * @brief Flag to define whether or not the grid has been prepared for integration
     */
    bool _grid_prepared;
    /**
     * @brief Flag to define whether or not the generation has been prepared using @a SetGen (very time-consuming operation, thus needs to be called once)
     */
    bool _generation_prepared;
    /**
     * @brief Total number of iterations for the current Vegas instance
     */
    int _mds;
    double _acc;
    double _alph;
    int _it;
    double _si;
    double _si2;
    double _swgt;
    double _schi;
    double _scalls;
    unsigned int _nd;
    unsigned int _ng;
    unsigned int _npg;
    double _calls;
    double _dxg;
    double _dv2g;
    double _xnd;
    unsigned int _ndm;
    double _xjac;
    int _now;
    double _vegas_result;
    double _vegas_abserr;
    /**
     * @brief The function which will be integrated by this Vegas instance
     * @param x_ The point at which this function is evaluated
     * @param ndim_ The number of degrees of freedom this function has
     * @param params_ A "_void_-ified" Parameters object to define the boundaries of the phase space (physics constraints)
     */
    double (*_f)(double* x_, size_t ndim_, void* params_);
};

#endif

