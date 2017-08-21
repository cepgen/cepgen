#include "StructureFunctions.h"

namespace CepGen
{
  std::ostream&
  operator<<( std::ostream& os, const StructureFunctions& sf )
  {
    switch ( sf ) {
      case Electron:            os << "electron"; break;
      case ElasticProton:       os << "elastic proton"; break;
      case SuriYennie:          os << "Suri-Yennie"; break;
      case SuriYennieLowQ2:     os << "Suri-Yennie;lowQ2"; break;
      case SzczurekUleshchenko: os << "Szczurek-Uleshchenko"; break;
      case FioreVal:            os << "Fiore;valence"; break;
      case FioreSea:            os << "Fiore;sea"; break;
      case Fiore:               os << "Fiore"; break;
    }
    return os;
  }
}
