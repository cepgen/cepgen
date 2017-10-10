#include "CepGen/IO/FortranInterface.h"

#include "CepGen/StructureFunctions/SuriYennie.h"
#include "CepGen/StructureFunctions/SzczurekUleshchenko.h"
#include "CepGen/StructureFunctions/FioreBrasse.h"
#include "CepGen/StructureFunctions/ChristyBosted.h"
#include "CepGen/StructureFunctions/BlockDurandHa.h"
#include "CepGen/StructureFunctions/ALLM.h"
#include "CepGen/StructureFunctions/Schaefer.h"

#include "MSTWGridHandler.h"

void
cepgen_structure_functions_( int& sfmode, double& q2, double& xbj, double& f2, double& fl )
{
  const double q2arg = q2, xbjarg = xbj;
  CepGen::StructureFunctions sf;
  switch ( ( SFmode )sfmode ) {
    case SuriYennie:          { const CepGen::SF::SuriYennie sy; sf = sy( q2arg, xbjarg ); } break;
    case SzczurekUleshchenko: { const CepGen::SF::SzczurekUleshchenko su; sf = su( q2arg, xbjarg ); } break;
    case BlockDurandHa:       { const CepGen::SF::BlockDurandHa bdh; sf = bdh( q2arg, xbjarg ); } break;
    case MSTWgrid: { sf = MSTW::GridHandler::get( "External/F2_Luxlike_fit/mstw_f2_scan_nnlo.txt" ).eval( q2, xbj ); } break;
    case ALLM91: { const CepGen::SF::ALLM allm( CepGen::SF::ALLM::Parameterisation::allm91() ); sf = allm( q2arg, xbjarg ); } break;
    case ALLM97: { const CepGen::SF::ALLM allm( CepGen::SF::ALLM::Parameterisation::allm97() ); sf = allm( q2arg, xbjarg ); } break;
    case GD07p:  { const CepGen::SF::ALLM allm( CepGen::SF::ALLM::Parameterisation::gd07p() ); sf = allm( q2arg, xbjarg ); } break;
    case GD11p:  { const CepGen::SF::ALLM allm( CepGen::SF::ALLM::Parameterisation::gd11p() ); sf = allm( q2arg, xbjarg ); } break;
    case Schaefer: { const CepGen::SF::Schaefer luxlike; sf = luxlike( q2arg, xbjarg ); } break;
    //--- resonances
    case FioreBrasse:   { const CepGen::SF::FioreBrasse fb; sf = fb( q2arg, xbjarg ); } break;
    case ChristyBosted: { const CepGen::SF::ChristyBosted cb; sf = cb( q2arg, xbjarg ); } break;
  };
  f2 = sf.F2;
  fl = sf.FL;
}

