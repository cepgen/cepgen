#ifndef _PHYSICS_H
#define _PHYSICS_H

#include "event.h"
#include "hadroniser.h"

extern "C"
{
  extern void grv95lo_(double&,double&,double&,double&,double&,double&,double&,double&);
}


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

Particles EPA(Particle* el_, Particle* pr_, int mode_, PhysicsBoundaries b_, double* q2_);

/**
 * Computes the proton structure function (F.W Brasse et al., DESY 76/11 (1976),
 *   http://dx.doi.org/10.1016/0550-3213(76)90231-5)
 * @cite Brasse1976413
 */
bool PSF(double,double,double*,double*,double*);

/**
 * @brief Vector meson particles and their decay mode
 */
typedef enum
{
  RHO_TO_PIPI = 113,
  OMEGA_TO_PIPI = 223,
  PHI_TO_KK = 333,
  PHI_TO_KLKS = 3332,
  JPSI_TO_LL = 444,
  PSIP_TO_LLX = 20443,
  UPS1S_TO_LL = 553,
  UPS2S_TO_LLX = 20553,
  UPS3S_TO_LLX = 30553,
  RHO1450_TO_PIPIRHO0 = 40113,
  PHI1680_TO_KKBAR = 10333
} VMDecay;

/**
 * Gets the branching ratio for a decay process, given its VMDecay identifier
 * @param[in] processId_ The identifier of the process
 * @return Branching ratio for the process
 */
double GetBRFromProcessId(Particle::ParticleCode vmId_);

/**
 * Decay a vector meson given its branching fractions
 * @param[in] part_ A Particle object containing all the physical and kinematic quantities for its decay
 * @param[in] had_ The default hadroniser object to use for the default decay
 * @return A vector of Particle objects
 */
Particles VMDecayer(Particle part_, Hadroniser *had_);

double ElasticFlux(double x_, double kt2_);
double InelasticFlux(double x_, double kt2_, double mx_);

#endif
