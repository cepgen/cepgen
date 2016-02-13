#ifndef Physics_h
#define Physics_h

#include "Event.h"
//#include "GenericHadroniser.h"

extern "C"
{
  //extern void grv95lo_(double&,double&,double&,double&,double&,double&,double&,double&);
  extern void grv95lo_(float&,float&,float&,float&,float&,float&,float&,float&);
}

class GenericHadroniser; // forward

/// List of physical constraints to apply on the phase space
class PhysicsBoundaries
{
 public:
  PhysicsBoundaries();
  ~PhysicsBoundaries();
  /// Minimal centre-of-mass energy for a \f$\gamma p\f$ system, in GeV.
  double wmin;
  /// Maximal centre-of-mass energy for a \f$\gamma p\f$ system, in GeV.
  double wmax;
  /// Minimal virtuality \f$Q^2\f$ of a photon in GeV\f${}^2\f$
  double q2min;
  /// Maximal virtuality \f$Q^2\f$ of a photon in GeV\f${}^2\f$
  double q2max;
  /// Minimal value of a generic scaling variable \f$\zeta\f$
  double zmin;
  /// Maximal value of a generic scaling variable \f$\zeta\f$
  double zmax;
};

Particles EPA(Particle* el_, Particle* pr_, int mode_, PhysicsBoundaries b_, double* q2_);

/**
 * Compute the proton structure function (F.W Brasse et al., DESY 76/11 (1976),
 *   http://dx.doi.org/10.1016/0550-3213(76)90231-5)
 * @cite Brasse1976413
 */
bool PSF(double,double,double*,double*,double*);

/// Vector meson particles and their decay mode
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
 * Get the branching ratio for a decay process, given its VMDecay identifier
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
//Particles VMDecayer(Particle part_, GenericHadroniser *had_);

double ElasticFlux(double x_, double kt2_);
double InelasticFlux(double x_, double kt2_, double mx_);

struct FormFactors {
  double FE;
  double FM;
};
FormFactors TrivialFormFactors();
FormFactors ElasticFormFactors(double q2, double mi2);
FormFactors SuriYennieFormFactors(double q2, double mi2, double mf2);

/**
 * Lorentz boost of a 4-vector (from CERNLIB)
 * @param pi_ Input 4-vector to boost
 * @param pf_ Output boosted 4-vector
 * @author L. Pape
 * @date 20 Aug 1975
 * @author Ian McLaren (mclareni), CERN/CN
 * @date 14 Feb 1996
 */
void Lorenb(double u_, const Particle::Momentum& ps_, double pi_[], double pf_[]);


#endif
