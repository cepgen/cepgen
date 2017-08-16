#ifndef CepGen_Physics_StructureFunctions_h
#define CepGen_Physics_StructureFunctions_h

#include <iostream>

namespace CepGen
{
  /// Proton structure function to be used in the outgoing state description
  enum StructureFunctions {
    Electron = 1,
    ElasticProton = 2,
    SuriYennie = 11,
    SuriYennieLowQ2 = 12,
    SzczurekUleshchenko = 15,
    FioreVal = 101,
    FioreSea = 102,
    Fiore = 103
  };
  /// Human-readable format of a structure function object
  std::ostream& operator<<( std::ostream& os, const StructureFunctions& sf );
}

#endif
