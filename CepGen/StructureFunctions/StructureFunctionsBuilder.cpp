#include "StructureFunctionsBuilder.h"

#include "CepGen/StructureFunctions/ALLM.h"
#include "CepGen/StructureFunctions/BlockDurandHa.h"
#include "CepGen/StructureFunctions/FioreBrasse.h"
#include "CepGen/StructureFunctions/ChristyBosted.h"
#include "CepGen/StructureFunctions/GenericLHAPDF.h"
#include "CepGen/StructureFunctions/SuriYennie.h"
#include "CepGen/StructureFunctions/SzczurekUleshchenko.h"

namespace CepGen
{
  StructureFunctions
  StructureFunctionsBuilder::get( const StructureFunctionsType& sf_type, double q2, double xbj )
  {
    switch ( sf_type ) {
      case StructureFunctionsType::SzczurekUleshchenko: { const SF::SzczurekUleshchenko su; return su( q2, xbj ); } break;
      case StructureFunctionsType::SuriYennie:          { const SF::SuriYennie sy; return sy( q2, xbj ); } break;
      case StructureFunctionsType::FioreBrasse:         { const SF::FioreBrasse fb; return fb( q2, xbj ); } break;
      case StructureFunctionsType::ChristyBosted:       { const SF::ChristyBosted cb; return cb( q2, xbj ); } break;
      case StructureFunctionsType::BlockDurandHa:       { const SF::BlockDurandHa bdh; return bdh( q2, xbj ); } break;
      case StructureFunctionsType::ALLM91:              { const SF::ALLM allm91( SF::ALLM::Parameterisation::allm91() ); return allm91( q2, xbj ); } break;
      case StructureFunctionsType::ALLM97:              { const SF::ALLM allm97( SF::ALLM::Parameterisation::allm97() ); return allm97( q2, xbj ); } break;
      case StructureFunctionsType::GD07p:               { const SF::ALLM gd07p( SF::ALLM::Parameterisation::gd07p() ); return gd07p( q2, xbj ); } break;
      case StructureFunctionsType::GD11p:               { const SF::ALLM gd11p( SF::ALLM::Parameterisation::gd11p() ); return gd11p( q2, xbj ); } break;
      default: return StructureFunctions(); //FIXME
    }
  }
}
