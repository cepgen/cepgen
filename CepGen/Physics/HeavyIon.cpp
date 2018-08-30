#include "CepGen/Physics/HeavyIon.h"

#include <sstream>

namespace CepGen
{
  HeavyIon
  HeavyIon::fromPDG( const PDG& pdg )
  {
    unsigned int ipdg = (unsigned int)pdg;
    return HeavyIon{ (unsigned short)( ipdg % 1000 ),
                     (Element)( ( ipdg / 1000 ) % 1000 ) };
  }

  HeavyIon::operator PDG() const
  {
    // (Pythia8 convention/10-1e10+1e6)
    return (PDG)( 1000000+1000*(unsigned short)Z+A );
  }

  std::ostream&
  operator<<( std::ostream& os, const Element& elem )
  {
    switch ( elem ) {
      case Element::H:  return os << "H";
      case Element::C:  return os << "C";
      case Element::O:  return os << "O";
      case Element::Cu: return os << "Cu";
      case Element::Xe: return os << "Xe";
      case Element::Au: return os << "Au";
      case Element::Pb: return os << "Pb";
    }
    return os;
  }

  std::ostream&
  operator<<( std::ostream& os, const HeavyIon& hi )
  {
    std::ostringstream oss; oss << hi.Z;
    if ( oss.str().empty() )
      return os << "HI{Z=" << hi.Z << ", A=" << hi.A << "}";
    return os << hi.A << oss.str();
  }
}
