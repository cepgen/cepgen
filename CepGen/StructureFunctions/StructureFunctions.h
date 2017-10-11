#ifndef CepGen_StructureFunctions_StructureFunctions_h
#define CepGen_StructureFunctions_StructureFunctions_h

#include <iostream>

namespace CepGen
{
  class StructureFunctions
  {
    public:
      StructureFunctions( double f2=0.0 ) : F2( f2 ), FL( 0.0 ), FM( 0.0 ) {}
      /// Proton structure function to be used in the outgoing state description
      /// \note Values correspond to the LPAIR legacy steering card values
      enum Type {
        Electron            = 1,
        ElasticProton       = 2,
        SuriYennie          = 11,
        SzczurekUleshchenko = 12,
        BlockDurandHa       = 13,
        FioreBrasse         = 101,
        ChristyBosted       = 102,
        ALLM91              = 201,
        ALLM97              = 202,
        GD07p               = 203,
        GD11p               = 204
      };

      double F2, F1;
      double FL, FM;
  };
  std::ostream& operator<<( std::ostream&, const StructureFunctions& );
  std::ostream& operator<<( std::ostream&, const StructureFunctions::Type& );
}

#endif
