#include "CepGen/Physics/HeavyIon.h"

namespace CepGen
{
  HeavyIon
  HeavyIon::fromPDG( const PDG& pdg )
  {
    unsigned int ipdg = (unsigned int)pdg;
    return HeavyIon{ (unsigned short)( ipdg % 1000 ), (unsigned short)( ( ipdg / 1000 ) % 1000 ) };
  }

  HeavyIon::operator PDG() const
  {
    // (Pythia8 convention/10-1e10+1e6)
    return (PDG)( 1e6+1e3*Z+A );
  }

  std::ostream&
  operator<<( std::ostream& os, const HeavyIon& hi )
  {
    switch ( hi.Z ) {
      case 6:   return os << hi.A << "C";
      case 8:   return os << hi.A << "O";
      case 29:  return os << hi.A << "Cu";
      case 54:  return os << hi.A << "Xe";
      case 79:  return os << hi.A << "Au";
      case 82:  return os << hi.A << "Pb";
      default:
        return os << "HI{Z=" << hi.Z << ", A=" << hi.A << "}";
    }
  }
}
