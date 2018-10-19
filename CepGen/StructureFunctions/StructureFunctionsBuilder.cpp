#include "CepGen/Core/ParametersList.h"

#include "CepGen/StructureFunctions/ALLM.h"
#include "CepGen/StructureFunctions/BlockDurandHa.h"
#include "CepGen/StructureFunctions/FioreBrasse.h"
#include "CepGen/StructureFunctions/ChristyBosted.h"
#include "CepGen/StructureFunctions/CLAS.h"
#include "CepGen/StructureFunctions/Partonic.h"
#include "CepGen/StructureFunctions/SuriYennie.h"
#include "CepGen/StructureFunctions/SzczurekUleshchenko.h"
#include "CepGen/StructureFunctions/Schaefer.h"
#include "CepGen/StructureFunctions/MSTWGrid.h"

namespace cepgen
{
  namespace strfun
  {
    std::shared_ptr<Parameterisation>
    Parameterisation::build( const ParametersList& params )
    {
      switch ( (Type)params.get<int>( "id" ) ) {
        case Type::Electron:
        case Type::ElasticProton:
        default:                        return std::make_shared<Parameterisation>();
        case Type::SzczurekUleshchenko: return std::make_shared<SzczurekUleshchenko>();
        case Type::SuriYennie:          return std::make_shared<SuriYennie>();
        case Type::FioreBrasse:         return std::make_shared<FioreBrasse>();
        case Type::ChristyBosted:       return std::make_shared<ChristyBosted>();
        case Type::CLAS:                return std::make_shared<CLAS>();
        case Type::BlockDurandHa:       return std::make_shared<BlockDurandHa>();
        case Type::ALLM91:              return std::make_shared<ALLM>( ParametersList().set<std::string>( "model", "ALLM91" ) );
        case Type::ALLM97:              return std::make_shared<ALLM>( ParametersList().set<std::string>( "model", "ALLM97" ) );
        case Type::GD07p:               return std::make_shared<ALLM>( ParametersList().set<std::string>( "model", "GD07p" ) );
        case Type::GD11p:               return std::make_shared<ALLM>( ParametersList().set<std::string>( "model", "GD11p" ) );
        case Type::Schaefer:            return std::make_shared<Schaefer>();
        case Type::Partonic:            return std::make_shared<Partonic>();
        //--- particular case for the MSTW grid as we are dealing
        //--- with a singleton ; hence, no deleter is needed!
        case Type::MSTWgrid:            return std::shared_ptr<mstw::Grid>( &mstw::Grid::get(), [=]( mstw::Grid* ){} );
      }
    }
  }
}
