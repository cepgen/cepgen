#include "StructureFunctionsBuilder.h"

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
  StructureFunctionsBuilder::get( const StructureFunctions::Type& sf_type, double q2, double xbj )
  {
    switch ( sf_type ) {
      case StructureFunctions::Electron:
      case StructureFunctions::ElasticProton:
        return StructureFunctions();
        break;
      case StructureFunctions::SzczurekUleshchenko: {
        const SF::SzczurekUleshchenko su;
        return su( q2, xbj );
      } break;
      case StructureFunctions::SuriYennie: {
        const SF::SuriYennie sy;
        return sy( q2, xbj );
      } break;
      case StructureFunctions::FioreBrasse: {
        const SF::FioreBrasse fb;
        return fb( q2, xbj );
      } break;
      case StructureFunctions::ChristyBosted: {
        const SF::ChristyBosted cb;
        return cb( q2, xbj );
      } break;
      case StructureFunctions::CLAS: {
        const SF::CLAS clas;
        return clas( q2, xbj );
      } break;
      case StructureFunctions::BlockDurandHa: {
        const SF::BlockDurandHa bdh;
        return bdh( q2, xbj );
      } break;
      case StructureFunctions::ALLM91: {
        const SF::ALLM allm91( SF::ALLM::Parameterisation::allm91() );
        return allm91( q2, xbj );
      } break;
      case StructureFunctions::ALLM97: {
        const SF::ALLM allm97( SF::ALLM::Parameterisation::allm97() );
        return allm97( q2, xbj );
      } break;
      case StructureFunctions::GD07p: {
        const SF::ALLM gd07p( SF::ALLM::Parameterisation::gd07p() );
        return gd07p( q2, xbj );
      } break;
      case StructureFunctions::GD11p: {
        const SF::ALLM gd11p( SF::ALLM::Parameterisation::gd11p() );
        return gd11p( q2, xbj );
      } break;
      case StructureFunctions::Schaefer: {
        const SF::Schaefer luxlike;
        return luxlike( q2, xbj );
      } break;
      case StructureFunctions::MSTWgrid: {
        return MSTW::GridHandler::get().eval( q2, xbj );
      } break;
    }
    return StructureFunctions(); //FIXME
  }
}
