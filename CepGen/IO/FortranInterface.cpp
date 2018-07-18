#include "CepGen/StructureFunctions/StructureFunctionsBuilder.h"
#include "CepGen/IO/MSTWGridHandler.h"

extern "C"
{
  void
  cepgen_structure_functions_( int& sfmode, double& q2, double& xbj, double& f2, double& fl )
  {
    using namespace CepGen;
    SF::Type sf_mode = (SF::Type)sfmode;
    if ( sf_mode == SF::Type::MSTWgrid ) {
      StructureFunctions sf = MSTW::GridHandler::get().eval( q2, xbj );
      f2 = sf.F2;
      fl = sf.FL;
      return;
    }
    StructureFunctions sf = StructureFunctionsBuilder::get( sf_mode )( q2, xbj );
    f2 = sf.F2;
    fl = sf.FL;
  }
}

