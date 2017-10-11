#include "CepGen/StructureFunctions/StructureFunctionsBuilder.h"

extern "C"
{
  void
  cepgen_structure_functions_( int& sfmode, double& q2, double& xbj, double& f2, double& fl )
  {
    const CepGen::StructureFunctions::Type mode = ( CepGen::StructureFunctions::Type )sfmode;
    const CepGen::StructureFunctions sf = CepGen::StructureFunctionsBuilder::get( mode, q2, xbj );
    f2 = sf.F2;
    fl = sf.FL;
  }
}

