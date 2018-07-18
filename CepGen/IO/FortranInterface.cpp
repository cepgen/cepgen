#include "CepGen/StructureFunctions/StructureFunctionsBuilder.h"
#include "CepGen/IO/MSTWGridHandler.h"
#include "CepGen/Core/Exception.h"

extern "C"
{
  void
  cepgen_structure_functions_( int& sfmode, double& q2, double& xbj, double& f2, double& fl )
  {
    using namespace CepGen;
    SF::Type sf_mode = (SF::Type)sfmode;

    CG_DEBUG( "cepgen_structure_functions" ) << sf_mode;

    if ( sf_mode == SF::Type::MSTWgrid ) {
      StructureFunctions sf = MSTW::GridHandler::get().eval( q2, xbj );
      f2 = sf.F2;
      fl = sf.FL;
      return;
    }
    StructureFunctions* sf = StructureFunctionsBuilder::get( sf_mode );
    StructureFunctions val = ( *sf )( q2, xbj );
    f2 = val.F2;
    fl = val.FL;
    delete sf;
  }
}

