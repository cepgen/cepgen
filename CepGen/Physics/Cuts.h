#ifndef CepGen_Physics_Cuts_h
#define CepGen_Physics_Cuts_h

#include "CepGen/Physics/Limits.h"
#include <vector>

namespace cepgen
{
  /// Constraints to be applied on the events kinematics
  struct Cuts
  {
    Limits pt_single;       ///< single particle transverse momentum
    Limits eta_single;      ///< single particle pseudo-rapidity
    Limits rapidity_single; ///< single particle rapidity
    Limits energy_single;   ///< single particle energy
    Limits mass_single;     ///< single particle mass
    Limits pt_sum;          ///< multiparticle system transverse momentum
    Limits eta_sum;         ///< multiparticle system pseudo-rapidity
    Limits energy_sum;      ///< multiparticle system energy
    Limits mass_sum;        ///< multiparticle system invariant mass
    Limits pt_diff;         ///< transverse momentum balance between the central particles
    Limits phi_pt_diff;     ///< azimuthal angles difference between the central particles
    Limits rapidity_diff;   ///< rapidity balance between the central particles
    Limits q2;              ///< parton virtuality
    Limits qt;              ///< parton transverse virtuality
    Limits phi_qt;          ///< parton azimuthal angle difference
    /// A collection of name -> limits
    std::vector<std::pair<std::string,Limits> > list() const;
  };
}

#endif
