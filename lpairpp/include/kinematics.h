#ifndef _KINEMATICS_H
#define _KINEMATICS_H

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>

#include "utils.h"


/**
 * @brief List of kinematic cuts to apply on the central and outgoing phase
 * space.
 */
class Kinematics
{
 public:
  Kinematics();
  ~Kinematics();
  /**
   * @brief Dumps all the parameters used in this process cross-section
   * computation / events generation
   */
  void Dump();
  /**
   * Type of kinematics to consider for the process. Can either be :
   *  * 0 for the electron-electron elastic case
   *  * 1 for the proton-proton elastic case
   *  * 2 for the proton-proton single-dissociative (or inelastic) case
   *  * 3 for the proton-proton double-dissociative case
   * @brief Type of kinematics to consider for the phase space
   */
  int kinematics;
  /**
   * @brief Sets of cuts to apply on the final phase space
   */
  int mode;
  /**
   * @brief Minimal transverse momentum of the single outgoing leptons
   */
  double ptmin;
  /**
   * @brief Maximal transverse momentum of the single outgoing leptons
   */
  double ptmax;
  /**
   * @brief Minimal energy of the central two-photons system
   */
  double emin;
  /**
   * @brief Maximal energy of the central two-photons system
   */
  double emax;
  /**
   * @brief Minimal polar (\f$\theta_\mathrm{min}\f$) angle of the outgoing
   * leptons, expressed in degrees
   */
  double thetamin;
  /**
   * @brief Maximal polar (\f$\theta_\mathrm{max}\f$) angle of the outgoing
   * leptons, expressed in degrees
   */
  double thetamax;
  /**
   * @brief Minimal mass (in GeV/c\f${}^\mathrm{2}\f$) of the outgoing proton
   * remnant(s)
   */
  double mxmin;
  /**
   * @brief Maximal mass (in GeV/c\f${}^\mathrm{2}\f$) of the outgoing proton
   * remnant(s)
   */
  double mxmax;
  /**
   * @brief The minimal value of \f$Q^2\f$
   */
  double q2min;
  /**
   * @brief The maximal value of \f$Q^2\f$
   */
  double q2max;
  /**
   * @brief The minimal \f$s\f$ on which the cross section is integrated
   */
  double wmin;
  /**
   * @brief The maximal \f$s\f$ on which the cross section is integrated. If
   * negative, the maximal energy available to the system (hence,
   * \f$s=(\sqrt{s})^{2}\f$) is provided.
   */
  double wmax;
};

#endif

