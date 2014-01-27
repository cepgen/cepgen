#ifndef _MCGen_H
#define _MCGen_H

#include <fstream>
#include <sstream>
#include <string>
#include <ctime>

#include "vegas.h"
#include "gamgam.h"
#include "parameters.h"

////////////////////////////////////////////////////////////////////////////////

/**
 * @image latex lpair_logo.pdf
 * @mainpage Principles
 * This Monte Carlo generator, based on the LPAIR code developed in the
 * early 1990s by J. Vermaseren *et al*@cite Vermaseren1983347, allows to compute
 * the cross-section and to generate events for the \f$\gamma\gamma\to\ell^{+}\ell^{-}\f$
 * process in high energy physics.
 * 
 * The main operation is the integration of the matrix element (given as a GamGam
 * object, subset of a Process object) performed by *Vegas*, an importance sampling
 * algorithm written in 1972 by G. P. Lepage@cite PeterLepage1978192. 
 *
 */

////////////////////////////////////////////////////////////////////////////////

/**
 * This object represents the core of this Monte Carlo generator, with its
 * allowance to generate the events (using the embedded Vegas object) and to
 * study the phase space in term of the variation of resulting cross section
 * while scanning the various parameters (point \f$\textbf{x}\f$ in the
 * DIM-dimensional phase space).
 *
 * The phase space is constrained using the InputParameters object given as an
 * argument to the constructor, and the differential cross-sections for each
 * value of the array \f$\textbf{x}\f$ are computed in the f-function defined
 * outside (but populated inside) this object.
 *
 * This f-function embeds a GamGam object which defines all the methods to
 * obtain this differential cross-section as well as the in- and outgoing
 * kinematics associated to each particle.
 *
 * @author Laurent Forthomme <laurent.forthomme@uclouvain.be>
 * @date February 2013
 * @brief Core of the Monte-Carlo generator
 *
 */
class MCGen {
 public:
  /**
   * Sets the number of dimensions on which to perform the integration, according to the set of input parameters given as an argument and propagated to the whole object
   * @brief Class constructor
   * @param[in] ip_ List of input parameters defining the phase space on which to perform the integration
   */
  MCGen(Parameters *ip_);
  ~MCGen();
  void Test();
  /**
   * Computes the cross-section for the run defined by this object. This returns the cross-section as well as the absolute error computed along.
   * @brief Compute the cross-section for the given process
   * @param[out] xsec_ The computed cross-section, in pb
   * @param[out] err_ The absolute integration error on the computed cross-section, in pb
   */
  void ComputeXsection(double* xsec_,double* err_);
  /**
   * Generates one single event given the phase space computed by Vegas in the integration step
   * @return A pointer to the Event object generated in this run
   */
  Event* GenerateOneEvent();
  /**
   * Launches the full events generation 
   * @deprecated This method is to be suppressed since the events generation can now be launched one event at a time using the @a GenerateOneEvent method
   */
  void LaunchGeneration();
  /**
   * @brief Returns the set of parameters used to setup the phase space to integrate
   * @return The Parameter object embedded in this class
   */
  Parameters GetParameters() { return *_par; }
  void AnalyzePhaseSpace(const std::string);
 private:
  /**
   * The GamGam object which allows to compute the outgoing particles' kinematics
   * as well as the cross-section for the given point in the phase space
   */
  GamGam *gg;
  /** @brief The Vegas integrator which will integrate the function */
  Vegas *veg;
  /** @brief Set of parameters to setup the phase space to integrate */
  Parameters *_par;
  double _xsec;
  double _xsec_error;
};

/**
 * The function to be integrated, which returns the value of the weight of an
 * event, including the matrix element of the process, all the kinematic
 * factors, and the cut restrictions. \f$x\f$ is an array of random numbers used
 * to select a random point inside the phase space.
 */
double f(double*,size_t,void*);

#endif

