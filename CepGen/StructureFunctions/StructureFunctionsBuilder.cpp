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
    Parameterisation::build( const Type& type, const ParametersList& params )
    {
      ParametersList pcopy = params;
      return build( pcopy.set<int>( "id", (int)type ) );
    }

    std::shared_ptr<Parameterisation>
    Parameterisation::build( const ParametersList& params )
    {
      ParametersList pcopy = params;
      switch ( (Type)params.get<int>( "id" ) ) {
        case Type::Electron:
        case Type::ElasticProton:
        default:                        return std::make_shared<Parameterisation>( params );
        case Type::SzczurekUleshchenko: return std::make_shared<SzczurekUleshchenko>( params );
        case Type::SuriYennie:          return std::make_shared<SuriYennie>( params );
        case Type::FioreBrasse:         return std::make_shared<FioreBrasse>( params );
        case Type::ChristyBosted:       return std::make_shared<ChristyBosted>( params );
        case Type::CLAS:                return std::make_shared<CLAS>( params );
        case Type::BlockDurandHa:       return std::make_shared<BlockDurandHa>( params );
        case Type::ALLM91:              return std::make_shared<ALLM>( pcopy.set<std::string>( "model", "ALLM91" ) );
        case Type::ALLM97:              return std::make_shared<ALLM>( pcopy.set<std::string>( "model", "ALLM97" ) );
        case Type::GD07p:               return std::make_shared<ALLM>( pcopy.set<std::string>( "model", "GD07p" ) );
        case Type::GD11p:               return std::make_shared<ALLM>( pcopy.set<std::string>( "model", "GD11p" ) );
        case Type::Schaefer:            return std::make_shared<Schaefer>( params );
        case Type::Partonic:            return std::make_shared<Partonic>( params );
        //--- particular case for the MSTW grid as we are dealing
        //--- with a singleton ; hence, no deleter is needed!
        case Type::MSTWgrid:            return std::shared_ptr<mstw::Grid>( mstw::Grid::get( params ).get(), [=]( mstw::Grid* ){} );
      }
    }
  }

  /// Human-readable format of a structure function type
  std::ostream&
  operator<<( std::ostream& os, const strfun::Type& sf )
  {
    switch ( sf ) {
      case strfun::Type::Invalid:             return os << "[INVALID]";
      case strfun::Type::Electron:            return os << "electron";
      case strfun::Type::ElasticProton:       return os << "elastic proton";
      case strfun::Type::SuriYennie:          return os << "Suri-Yennie";
      case strfun::Type::SzczurekUleshchenko: return os << "Szczurek-Uleshchenko";
      case strfun::Type::FioreBrasse:         return os << "Fiore-Brasse";
      case strfun::Type::ChristyBosted:       return os << "Christy-Bosted";
      case strfun::Type::CLAS:                return os << "CLAS";
      case strfun::Type::BlockDurandHa:       return os << "BDH";
      case strfun::Type::ALLM91:              return os << "ALLM91";
      case strfun::Type::ALLM97:              return os << "ALLM97";
      case strfun::Type::GD07p:               return os << "GD07p";
      case strfun::Type::GD11p:               return os << "GD11p";
      case strfun::Type::Schaefer:            return os << "LUXlike";
      case strfun::Type::MSTWgrid:            return os << "MSTW (grid)";
      case strfun::Type::Partonic:            return os << "Partonic";
    }
    return os;
  }
}
