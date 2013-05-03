#ifndef _MCGen_H
#define _MCGen_H

#include <fstream>
#include <sstream>
#include <string>
#include "vegas.h"
#include "utils.h"
#include "gamgam.h"
#include "gnuplot.h"

/** @brief Core Monte-Carlo generator */
class MCGen {
 public:
  /**
   * Sets the number of dimensions on which to perform the integration, according
   *  to the set of input parameters given as an argument and propagated to the
   *  whole object
   * @brief Class constructor
   * @param ip_ List of input parameters defining the phase space
   *  on which to perform the integration
   */
  MCGen(InputParameters);
  ~MCGen();
  void LaunchGen(int);
  /** 
   * @brief Returns the set of parameters used to setup the phase space to
   *   integrate
   * @return The InputParameter object embedded in this class
   */
  InputParameters GetInputParameters() { return _ip; }
  void AnalyzePhaseSpace(std::string);
 private:
  /**
   * @brief The GamGam object computing the kinematics and the cross-section for the
   * given point in the phase space
   */
  GamGam *gg;
  /** @brief The Vegas integrator which will integrate the function */
  Vegas *veg;
  /** @brief Number of dimensions on which to perform the integration */
  int _ndim;
  int _inp1pdg, _inp2pdg, _outp1pdg, _outp2pdg;
  /** @brief Set of parameters to setup the phase space to integrate */
  InputParameters _ip;  
};

/**
 * The function to be integrated, which returns the value of the weight of an
 * event, including the matrix element of the process, all the kinematic
 * factors, and the cut restrictions. \f$x\f$ is an array of random numbers used
 * to select a random point inside the phase space.
 */
double f(double*,size_t,void*);

#endif

