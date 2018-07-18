#include "CepGen/StructureFunctions/StructureFunctionsBuilder.h"
#include "CepGen/StructureFunctions/StructureFunctions.h"
#include "CepGen/Core/Exception.h"

#ifdef __cplusplus
extern "C" {
#endif
  void
  cepgen_structure_functions_( int& sfmode, double& q2, double& xbj, double& f2, double& fl )
  {
    using namespace CepGen;
    SF::Type sf_mode = (SF::Type)sfmode;

    CG_DEBUG( "cepgen_structure_functions" ) << sf_mode;

    StructureFunctions& val = ( *StructureFunctionsBuilder::get( sf_mode ) )( q2, xbj );
    f2 = val.F2;
    fl = val.FL;
  }
#ifdef __cplusplus
}
#endif

