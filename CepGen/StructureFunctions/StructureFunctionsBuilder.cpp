#include "CepGen/StructureFunctions/StructureFunctionsBuilder.h"

#include "CepGen/StructureFunctions/ALLM.h"
#include "CepGen/StructureFunctions/BlockDurandHa.h"
#include "CepGen/StructureFunctions/FioreBrasse.h"
#include "CepGen/StructureFunctions/ChristyBosted.h"
#include "CepGen/StructureFunctions/CLAS.h"
#include "CepGen/StructureFunctions/GenericLHAPDF.h"
#include "CepGen/StructureFunctions/SuriYennie.h"
#include "CepGen/StructureFunctions/SzczurekUleshchenko.h"
#include "CepGen/StructureFunctions/Schaefer.h"
#include "CepGen/IO/MSTWGridHandler.h"

namespace CepGen
{
  StructureFunctions
  StructureFunctionsBuilder::get( const SF::Type& sf_type )
  {
    std::cout << sf_type << std::endl;
    switch ( sf_type ) {
      case SF::Type::Electron:
      case SF::Type::MSTWgrid: //FIXME
      default:
      case SF::Type::ElasticProton:       return StructureFunctions();
      case SF::Type::SzczurekUleshchenko: return SF::SzczurekUleshchenko();
      case SF::Type::SuriYennie:          return SF::SuriYennie();
      case SF::Type::FioreBrasse:         return SF::FioreBrasse();
      case SF::Type::ChristyBosted:       return SF::ChristyBosted();
      case SF::Type::CLAS:                return SF::CLAS();
      case SF::Type::BlockDurandHa:       return SF::BlockDurandHa();
      case SF::Type::ALLM91:              return SF::ALLM( SF::ALLM::Parameterisation::allm91() );
      case SF::Type::ALLM97:              return SF::ALLM( SF::ALLM::Parameterisation::allm97() );
      case SF::Type::GD07p:               return SF::ALLM( SF::ALLM::Parameterisation::gd07p() );
      case SF::Type::GD11p:               return SF::ALLM( SF::ALLM::Parameterisation::gd11p() );
      case SF::Type::Schaefer:            return SF::Schaefer();
    }
  }
}
