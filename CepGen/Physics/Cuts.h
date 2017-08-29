#ifndef CepGen_Physics_Cuts_h
#define CepGen_Physics_Cuts_h

namespace CepGen
{
  /// Constraints to be applied on the events kinematics
  namespace Cuts
  {
    /// Cuts on the central particles (e.g. dilepton system for LPAIR)
    enum Central
    {
      pt_single,     ///< single particle transverse momentum
      eta_single,    ///< single particle pseudo-rapidity
      energy_single, ///< single particle energy
      mass_sum,      ///< central system invariant mass
      pt_sum,        ///< central system transverse momentum
      eta_sum,       ///< central system pseudo-rapidity
      energy_sum,    ///< central system energy
      pt_diff,       ///< transverse momentum balance between the central particles
      dely           ///< rapidity balance between the central particles
    };
    inline std::ostream& operator<<( std::ostream& os, const Central& is ) {
      switch ( is ) {
        case pt_single: return os << "Single central particle pT";
        case eta_single: return os << "Single central particle eta";
        case energy_single: return os << "Single central particle energy";
        case mass_sum: return os << "Central system mass";
        case pt_sum: return os << "Central system pT";
        case eta_sum: return os << "Central system eta";
        case energy_sum: return os << "Central system energy";
        case pt_diff: return os << "Central system d(pT)";
        case dely: return os << "Central system d(Y)";
      }
      return os;
    }

    /// Cuts on the beam particles remnants
    enum Remnants
    {
      mass ///< outgoing beam particles remnants mass
    };
    inline std::ostream& operator<<( std::ostream& os, const Remnants& is ) {
      switch ( is ) {
        case mass: return os << "Remnants mass";
      }
      return os;
    }

    /// Cuts on the initial state dynamics
    enum InitialState
    {
      q2, ///< parton virtuality
      qt, ///< parton transverse virtuality
      w   ///< two-parton squared momentum
    };
    inline std::ostream& operator<<( std::ostream& os, const InitialState& is ) {
      switch ( is ) {
        case q2: return os << "Virtuality range";
        case qt: return os << "Transverse virtuality range";
        case w: return os << "W";
      }
      return os;
    }
  }
}

#endif
