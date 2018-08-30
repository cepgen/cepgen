#include "CepGen/Physics/HeavyIon.h"

#include <sstream>

namespace CepGen
{
  HeavyIon::HeavyIon( const PDG& pdg ) :
    A( (unsigned int)pdg % 1000 ),
    Z( (Element)( (unsigned int)pdg/1000000 == 0 ? 0
                : ( (unsigned int)pdg/1000 ) % 1000 ) )
  {}

  HeavyIon::operator PDG() const
  {
    // Pythia8 convention/10-1e10+1e6
    return (PDG)( 1000000+1000*(unsigned short)Z+A );
  }

  std::ostream&
  operator<<( std::ostream& os, const HeavyIon& hi )
  {
    std::ostringstream oss; oss << hi.Z;
    if ( oss.str().empty() )
      return os << "HI{Z=" << (unsigned short)hi.Z << ", A=" << hi.A << "}";
    return os << hi.A << oss.str();
  }

  std::ostream&
  operator<<( std::ostream& os, const Element& elem )
  {
    switch ( elem ) {
      case Element::invalid:
        return os;
      case Element::H:  return os << "H";
      case Element::C:  return os << "C";
      case Element::O:  return os << "O";
      case Element::Al: return os << "Al";
      case Element::Cu: return os << "Cu";
      case Element::Xe: return os << "Xe";
      case Element::Au: return os << "Au";
      case Element::Pb: return os << "Pb";
      case Element::U:  return os << "U";
    }
    return os;
  }
}
