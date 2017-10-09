#include "CepGen/IO/FortranInterface.h"

#include "CepGen/StructureFunctions/SuriYennie.h"
#include "CepGen/StructureFunctions/SzczurekUleshchenko.h"
#include "CepGen/StructureFunctions/FioreBrasse.h"
#include "CepGen/StructureFunctions/ChristyBosted.h"
#include "CepGen/StructureFunctions/BlockDurandHa.h"
#include "CepGen/StructureFunctions/ALLM.h"
#include "CepGen/StructureFunctions/Schaefer.h"

void
cepgen_structure_functions_( int& sfmode, double& q2, double& xbj, double& f2, double& fl )
{
  const double q2arg = q2, xbjarg = xbj;
  switch ( ( SFmode )sfmode )
  {
    case SuriYennie: {
      const CepGen::SF::SuriYennie sy;
      const CepGen::StructureFunctions sf = sy( q2arg, xbjarg );
      f2 = sf.F2;
      fl = sf.FL;
    } break;
    case SzczurekUleshchenko: {
      const CepGen::SF::SzczurekUleshchenko su;
      const CepGen::StructureFunctions sf = su( q2arg, xbjarg );
      f2 = sf.F2;
      fl = 0.;
    } break;
    case BlockDurandHa: {
      const CepGen::SF::BlockDurandHa bdh;
      const CepGen::StructureFunctions sf = bdh( q2arg, xbjarg );
      f2 = sf.F2;
      fl = 0.;
    } break;
    case ALLM91: {
      const CepGen::SF::ALLM allm( CepGen::SF::ALLM::Parameterisation::allm91() );
      const CepGen::StructureFunctions sf = allm( q2arg, xbjarg );
      f2 = sf.F2;
      fl = 0.;
    } break;
    case ALLM97: {
      const CepGen::SF::ALLM allm( CepGen::SF::ALLM::Parameterisation::allm97() );
      const CepGen::StructureFunctions sf = allm( q2arg, xbjarg );
      f2 = sf.F2;
      fl = 0.;
    } break;
    case GD07p: {
      const CepGen::SF::ALLM allm( CepGen::SF::ALLM::Parameterisation::gd07p() );
      const CepGen::StructureFunctions sf = allm( q2arg, xbjarg );
      f2 = sf.F2;
      fl = 0.;
    } break;
    case GD11p: {
      const CepGen::SF::ALLM allm( CepGen::SF::ALLM::Parameterisation::gd11p() );
      const CepGen::StructureFunctions sf = allm( q2arg, xbjarg );
      f2 = sf.F2;
      fl = 0.;
    } break;
    case Schaefer: {
      const CepGen::SF::Schaefer luxlike;
      const CepGen::StructureFunctions sf = luxlike( q2arg, xbjarg );
      f2 = sf.F2;
      fl = sf.FL;
    } break;
    //--- resonances
    case FioreBrasse: {
      const CepGen::SF::FioreBrasse fb;
      const CepGen::StructureFunctions sf = fb( q2arg, xbjarg );
      f2 = sf.F2;
      fl = 0.;
    } break;
    case ChristyBosted: {
      const CepGen::SF::ChristyBosted cb;
      const CepGen::StructureFunctions sf = cb( q2arg, xbjarg );
      f2 = sf.F2;
      fl = 0.;
    } break;
  };
}

