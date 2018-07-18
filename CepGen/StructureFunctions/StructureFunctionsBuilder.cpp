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
#include "CepGen/StructureFunctions/MSTWGrid.h"

namespace CepGen
{
  StructureFunctions*
  StructureFunctionsBuilder::get( const SF::Type& sf_type )
  {
    switch ( sf_type ) {
      case SF::Type::Electron:
      case SF::Type::ElasticProton:
      default:                            return new StructureFunctions();
      case SF::Type::SzczurekUleshchenko: return new SF::SzczurekUleshchenko();
      case SF::Type::SuriYennie:          return new SF::SuriYennie();
      case SF::Type::FioreBrasse:         return new SF::FioreBrasse();
      case SF::Type::ChristyBosted:       return new SF::ChristyBosted();
      case SF::Type::CLAS:                return new SF::CLAS();
      case SF::Type::BlockDurandHa:       return new SF::BlockDurandHa();
      case SF::Type::ALLM91:              return new SF::ALLM( SF::ALLM::Parameterisation::allm91() );
      case SF::Type::ALLM97:              return new SF::ALLM( SF::ALLM::Parameterisation::allm97() );
      case SF::Type::GD07p:               return new SF::ALLM( SF::ALLM::Parameterisation::gd07p() );
      case SF::Type::GD11p:               return new SF::ALLM( SF::ALLM::Parameterisation::gd11p() );
      case SF::Type::Schaefer:            return new SF::Schaefer();
      case SF::Type::GenericLHAPDF:       return new SF::GenericLHAPDF();
      case SF::Type::MSTWgrid:            return &MSTW::Grid::get();
    }
  }
}
