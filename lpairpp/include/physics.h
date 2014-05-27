#ifndef _PHYSICS_H
#define _PHYSICS_H

#include "event.h"

/**
 * List of physical constraints to apply on the phase space
 */
class PhysicsBoundaries
{
 public:
  PhysicsBoundaries();
  ~PhysicsBoundaries();
  /**
   * @brief Minimal centre-of-mass energy for a \f$\gamma p\f$ system, in GeV.
   */
  double wmin;
  /**
   * @brief Maximal centre-of-mass energy for a \f$\gamma p\f$ system, in GeV.
   */
  double wmax;
  /**
   * @brief Minimal virtuality \f$Q^2\f$ of a photon in GeV\f${}^2\f$
   */
  double q2min;
  /**
   * @brief Maximal virtuality \f$Q^2\f$ of a photon in GeV\f${}^2\f$
   */
  double q2max;
  /**
   * @brief Minimal value of a generic scaling variable \f$\zeta\f$
   */
  double zmin;
  /**
   * @brief Maximal value of a generic scaling variable \f$\zeta\f$
   */
  double zmax;
};

Particles EPA(Particle el_, Particle pr_, int mode_, PhysicsBoundaries b_);

/**
 * Computes the proton structure function (F.W Brasse et al., DESY 76/11 (1976),
 *   http://dx.doi.org/10.1016/0550-3213(76)90231-5)
 * @cite Brasse1976413
 */
bool PSF(double,double,double*,double*,double*);

#endif
