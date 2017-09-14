#include "StructureFunctions.h"

namespace CepGen
{
  std::ostream&
  operator<<( std::ostream& os, const StructureFunctionsType& sf )
  {
    switch ( sf ) {
      case Electron:            return os << "electron";
      case ElasticProton:       return os << "elastic proton";
      case SuriYennie:          return os << "Suri-Yennie";
      case SuriYennieLowQ2:     return os << "Suri-Yennie;lowQ2";
      case SzczurekUleshchenko: return os << "Szczurek-Uleshchenko";
      case FioreVal:            return os << "Fiore;valence";
      case FioreSea:            return os << "Fiore;sea";
      case Fiore:               return os << "Fiore";
      case ALLM91:              return os << "ALLM;91";
      case ALLM97:              return os << "ALLM;97";
      case ALLM_HHT:            return os << "ALLM;HHT";
      case ALLM_HHT_FT:         return os << "ALLM;HHT-FT";
      default: return os;
    }
  }

  std::ostream&
  operator<<( std::ostream& os, const StructureFunctions& sf )
  {
    return os << ", F2 = " << sf.F2;
  }
}
