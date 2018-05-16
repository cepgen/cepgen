#ifndef CepGen_Physics_Cuts_h
#define CepGen_Physics_Cuts_h

namespace CepGen
{
  /// Constraints to be applied on the events kinematics
  enum class Cuts : unsigned short
  {
    pt_single,       ///< single particle transverse momentum
    eta_single,      ///< single particle pseudo-rapidity
    rapidity_single, ///< single particle rapidity
    energy_single,   ///< single particle energy
    mass_single,     ///< single particle mass
    pt_sum,          ///< multiparticle system transverse momentum
    eta_sum,         ///< multiparticle system pseudo-rapidity
    energy_sum,      ///< multiparticle system energy
    mass_sum,        ///< multiparticle system invariant mass
    pt_diff,         ///< transverse momentum balance between the central particles
    phi_pt_diff,     ///< azimuthal angles difference between the central particles
    rapidity_diff,   ///< rapidity balance between the central particles
    q2,              ///< parton virtuality
    qt,              ///< parton transverse virtuality
    phi_qt,          ///< parton azimuthal angle difference
    w                ///< two-parton squared momentum
  };
  inline std::ostream& operator<<( std::ostream& os, const Cuts& is ) {
    switch ( is ) {
      case Cuts::pt_single: return os << "Single central pt  (GeV/c)";
      case Cuts::eta_single: return os << "Single central eta";
      case Cuts::rapidity_single: return os << "Single central rapidity";
      case Cuts::energy_single: return os << "Single central energy (GeV)";
      case Cuts::mass_single: return os << "Single particle mass (GeV/c²)";
      case Cuts::pt_sum: return os << "Central system pt (GeV/c)";
      case Cuts::eta_sum: return os << "Central system eta";
      case Cuts::energy_sum: return os << "Central system energy";
      case Cuts::mass_sum: return os << "Central system mass";
      case Cuts::pt_diff: return os << "Central system Δpt (GeV/c)";
      case Cuts::phi_pt_diff: return os << "Central system Δɸ";
      case Cuts::rapidity_diff: return os << "Central system ΔY";
      case Cuts::q2: return os << "Virtuality range (GeV²)";
      case Cuts::qt: return os << "Transverse virtuality range (GeV)";
      case Cuts::phi_qt: return os << "Partons Δɸ range";
      case Cuts::w: return os << "W (GeV²)";
    }
    return os;
  }
}

#endif
