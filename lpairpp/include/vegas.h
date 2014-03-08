#ifndef _VEGAS_H
#define _VEGAS_H

#include <fstream>
#include <cstdio> // remove (DEBUG)

#include "parameters.h"

#define MAX_ND 50

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
     * @date September 1976
     * @date Reviewed in Apr 1978
     * @date FTN5 version 21 Aug 1984
     * @date This C++ implementation is from 12 Dec 2013
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
    double Treat(double* x_,Parameters* ip_,bool storedbg_=false);
    inline double Treat(double* x_) { return this->Treat(x_,(Parameters*)this->_ip); }
    inline double Treat(double* x_,bool storedbg_) { return this->Treat(x_,(Parameters*)this->_ip,storedbg_); }
    inline double F(double* x_) { return this->_f(x_, this->_ndim, (void*)this->_ip); }
    inline double F(double* x_,Parameters* ip_) { return this->_f(x_, this->_ndim, (void*)ip_); }
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
     * Sets all the generation mode variables and align them to the integration 
     * grid set while computing the cross-section
     * @brief Prepare the class for events generation
     */
    void SetGen();
    void DumpGrid();
    /** @brief The number of dimensions on which to integrate the function */
    const size_t _ndim;
    double (*_f)(double* x_, size_t ndim_, void* params_);
    unsigned int _ndo;
    int _nTreatCalls;
    int _nTreat;
    double _rTreat;
    double _mbin;
    double _ffmax;
    int *_n;
    int *_nm;
    double *_fmax;
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
    bool _force_correction;
    double *_xi[MAX_ND];
    double *_d[MAX_ND];
    double *_di[MAX_ND];
    Parameters *_ip;
    bool _grid_prepared;
};

#endif

