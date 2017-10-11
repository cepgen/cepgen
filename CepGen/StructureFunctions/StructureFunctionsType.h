#ifndef CepGen_StructureFunctions_StructureFunctionsType_h
#define CepGen_StructureFunctions_StructureFunctionsType_h

namespace CepGen
{
  /// Proton structure function to be used in the outgoing state description
  /// \note Values correspond to the LPAIR legacy steering card values
  enum StructureFunctionsType {
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
}

#endif
