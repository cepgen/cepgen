#include "CepGen/StructureFunctions/SzczurekUleshchenko.h"
#include "CepGen/StructureFunctions/SuriYennie.h"
#include "CepGen/StructureFunctions/FioreBrasse.h"
#include "CepGen/StructureFunctions/ChristyBosted.h"
#include "CepGen/StructureFunctions/CLAS.h"
#include "CepGen/StructureFunctions/BlockDurandHa.h"
#include "CepGen/StructureFunctions/ALLM.h"
#include "CepGen/StructureFunctions/Schaefer.h"
#include "CepGen/StructureFunctions/GenericLHAPDF.h"
#include "CepGen/IO/MSTWGridHandler.h"

extern "C"
{
  void
  cepgen_structure_functions_( int& sfmode, double& q2, double& xbj, double& f2, double& fl )
  {
    using namespace CepGen;
    StructureFunctions sf;
    switch ( (SF::Type)sfmode ) {
      case SF::Type::Electron:
      case SF::Type::ElasticProton:
      case SF::Type::Invalid: {
        f2 = fl = 0.;
        return;
      }
      case SF::Type::SzczurekUleshchenko: sf = SF::SzczurekUleshchenko(); break;
      case SF::Type::SuriYennie:          sf = SF::SuriYennie(); break;
      case SF::Type::FioreBrasse:         sf = SF::FioreBrasse(); break;
      case SF::Type::ChristyBosted:       sf = SF::ChristyBosted(); break;
      case SF::Type::CLAS:                sf = SF::CLAS(); break;
      case SF::Type::BlockDurandHa:       sf = SF::BlockDurandHa(); break;
      case SF::Type::ALLM91:              sf = SF::ALLM( SF::ALLM::Parameterisation::allm91() ); break;
      case SF::Type::ALLM97:              sf = SF::ALLM( SF::ALLM::Parameterisation::allm97() ); break;
      case SF::Type::GD07p:               sf = SF::ALLM( SF::ALLM::Parameterisation::gd07p() ); break;
      case SF::Type::GD11p:               sf = SF::ALLM( SF::ALLM::Parameterisation::gd11p() ); break;
      case SF::Type::Schaefer:            sf = SF::Schaefer(); break;
      case SF::Type::GenericLHAPDF:       sf = SF::GenericLHAPDF(); break;
      case SF::Type::MSTWgrid: {
        sf = MSTW::GridHandler::get().eval( q2, xbj );
        f2 = sf.F2;
        fl = sf.FL;
        return;
      }
    }
    sf = sf( q2, xbj );
    f2 = sf.F2;
    fl = sf.FL;
  }
}

