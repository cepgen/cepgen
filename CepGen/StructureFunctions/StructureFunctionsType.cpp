#include "StructureFunctionsType.h"

#include <iostream>

namespace CepGen
{
  /// Human-readable format of a structure function type
  std::ostream&
  operator<<( std::ostream& os, const StructureFunctionsType& sf )
  {
    switch ( sf ) {
      case Electron:            return os << "electron";
      case ElasticProton:       return os << "elastic proton";
      case SuriYennie:          return os << "Suri-Yennie";
      case SzczurekUleshchenko: return os << "Szczurek-Uleshchenko";
      case FioreBrasse:         return os << "Fiore-Brasse";
      case ChristyBosted:       return os << "Christy-Bosted";
      case BlockDurandHa:       return os << "BDH";
      case ALLM91:              return os << "ALLM;91";
      case ALLM97:              return os << "ALLM;97";
      case GD07p:               return os << "ALLM;GD07p";
      case GD11p:               return os << "ALLM;GD11p";
    }
    return os;
  }
}
