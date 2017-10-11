#ifndef CepGen_StructureFunctions_StructureFunctions_h
#define CepGen_StructureFunctions_StructureFunctions_h

#include "CepGen/Event/Particle.h"
#include <ostream>
#include <vector>
#include <cmath>

namespace CepGen
{
  /// Proton structure function to be used in the outgoing state description
  /// \note Values correspond to the LPAIR legacy steering card values
  enum StructureFunctionsType {
    Electron = 1,
    ElasticProton = 2,
    SuriYennie = 11,
    SuriYennieLowQ2 = 12,
    SzczurekUleshchenko = 15,
    FioreBrasse = 101,
    ALLM91, ALLM97, ALLM_HHT, ALLM_HHT_FT, ALLM_GD07p, ALLM_GD11p
  };
  /// Human-readable format of a structure function type
  std::ostream& operator<<( std::ostream& os, const StructureFunctionsType& sf );

  class StructureFunctions
  {
    public:
      StructureFunctions( double f2=0.0 ) : F2( f2 ), FL( 0.0 ), FM( 0.0 ) {}

      double F2;
      double F1;
      double FL, FM;
  };
  /// Human-readable format of a structure function object
  std::ostream& operator<<( std::ostream& os, const StructureFunctions& sf );
}

#endif
