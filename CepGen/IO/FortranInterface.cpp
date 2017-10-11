#include "CepGen/IO/FortranInterface.h"

#include "CepGen/StructureFunctions/StructureFunctionsBuilder.h"

void
cepgen_structure_functions_( int& sfmode, double& q2, double& xbj, double& f2, double& fl )
{
  const CepGen::StructureFunctions sf = CepGen::StructureFunctionsBuilder::get( ( CepGen::StructureFunctions::Type )sfmode, q2, xbj );
  f2 = sf.F2;
  fl = sf.FL;
}

