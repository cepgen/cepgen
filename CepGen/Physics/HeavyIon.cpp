#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/PDG.h"

#include "CepGen/Core/Exception.h"

#include <sstream>

namespace cepgen
{
  HeavyIon::HeavyIon( pdgid_t pdg ) :
    Z( (Element)( pdg/1000000 == 0 ? 0
                : ( pdg/1000 ) % 1000 ) ),
    A( ( Z != Element::invalid ) ? pdg % 1000 : 0 )
  {}

  HeavyIon::operator pdgid_t() const
  {
    // Pythia8 convention/10-1e10+1e6
    return (pdgid_t)( 1000000+1000*(unsigned short)Z+A );
  }

  HeavyIon::operator bool() const
  {
    return Z != Element::invalid; // skip the proton
  }

  double
  HeavyIon::mass( const HeavyIon& hi )
  {
    if ( !hi )
      throw CG_FATAL( "mass" ) << "Invalid heavy ion: " << hi << "!";
    return (short)hi.Z*PDG::get().mass( PDG::proton );
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
