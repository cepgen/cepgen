#include "CepGen/StructureFunctions/StructureFunctionsBuilder.h"

#include "CepGen/StructureFunctions/ALLM.h"
#include "CepGen/StructureFunctions/BlockDurandHa.h"
#include "CepGen/StructureFunctions/FioreBrasse.h"
#include "CepGen/StructureFunctions/ChristyBosted.h"
#include "CepGen/StructureFunctions/CLAS.h"
#include "CepGen/StructureFunctions/LHAPDF.h"
#include "CepGen/StructureFunctions/SuriYennie.h"
#include "CepGen/StructureFunctions/SzczurekUleshchenko.h"
#include "CepGen/StructureFunctions/Schaefer.h"
#include "CepGen/StructureFunctions/MSTWGrid.h"

namespace CepGen
{
  std::shared_ptr<StructureFunctions>
  StructureFunctionsBuilder::get( const SF::Type& sf_type )
  {
    switch ( sf_type ) {
      case SF::Type::Electron:
      case SF::Type::ElasticProton:
      default:                            return std::make_shared<StructureFunctions>();
      case SF::Type::SzczurekUleshchenko: return std::make_shared<SF::SzczurekUleshchenko>();
      case SF::Type::SuriYennie:          return std::make_shared<SF::SuriYennie>();
      case SF::Type::FioreBrasse:         return std::make_shared<SF::FioreBrasse>();
      case SF::Type::ChristyBosted:       return std::make_shared<SF::ChristyBosted>();
      case SF::Type::CLAS:                return std::make_shared<SF::CLAS>();
      case SF::Type::BlockDurandHa:       return std::make_shared<SF::BlockDurandHa>();
      case SF::Type::ALLM91:              return std::make_shared<SF::ALLM>( SF::ALLM::Parameterisation::allm91() );
      case SF::Type::ALLM97:              return std::make_shared<SF::ALLM>( SF::ALLM::Parameterisation::allm97() );
      case SF::Type::GD07p:               return std::make_shared<SF::ALLM>( SF::ALLM::Parameterisation::gd07p() );
      case SF::Type::GD11p:               return std::make_shared<SF::ALLM>( SF::ALLM::Parameterisation::gd11p() );
      case SF::Type::Schaefer:            return std::make_shared<SF::Schaefer>();
      case SF::Type::LHAPDF:              return std::make_shared<SF::LHAPDF>();
      //--- particular case for the MSTW grid as we are dealing
      //--- with a singleton ; hence, no deleter is needed!
      case SF::Type::MSTWgrid:            return std::shared_ptr<mstw::Grid>( &mstw::Grid::get(), [=]( mstw::Grid* ){} );
    }
  }
}
