#ifndef CepGen_Physics_Cuts_h
#define CepGen_Physics_Cuts_h

namespace CepGen
{
  namespace Cuts
  {
    enum Central
    {
      pt_single, eta_single, energy_single,
      mass_sum, pt_sum, eta_sum, energy_sum,
      pt_diff, dely
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


    enum Remnants
    {
      mass
    };
    inline std::ostream& operator<<( std::ostream& os, const Remnants& is ) {
      switch ( is ) {
        case mass: return os << "Remnants mass";
      }
      return os;
    }

    enum InitialState
    {
      q2, qt, w
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
