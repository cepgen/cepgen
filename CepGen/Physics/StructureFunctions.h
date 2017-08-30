#ifndef CepGen_Physics_StructureFunctions_h
#define CepGen_Physics_StructureFunctions_h

#include "Particle.h"

extern "C"
{
  extern void grv95lo_( float&, float&, float&, float&, float&, float&, float&, float& );
}

namespace CepGen
{
  /// Proton structure function to be used in the outgoing state description
  enum StructureFunctions {
    Electron = 1,
    ElasticProton = 2,
    SuriYennie = 11,
    SuriYennieLowQ2 = 12,
    SzczurekUleshchenko = 15,
    FioreVal = 101,
    FioreSea = 102,
    Fiore = 103
  };
  /// Human-readable format of a structure function object
  std::ostream& operator<<( std::ostream& os, const StructureFunctions& sf );

  /**
   * Compute the proton structure function (F.W Brasse et al., DESY 76/11 (1976),
   *   http://dx.doi.org/10.1016/0550-3213(76)90231-5)
   * \param[in] q2 Squared 4-momentum transfer
   * \param[in] mx2 Squared mass of the proton remnant
   * \param[out] sigma_t ...
   * \param[out] w1 First proton structure function: \f$\mathcal W_1\f$
   * \param[out] w2 Second proton structure function: \f$\mathcal W_2\f$
   * \cite Brasse1976413
   */
  bool PSF( double q2, double mx2, double& sigma_t, double& w1, double& w2 );

}

#endif
