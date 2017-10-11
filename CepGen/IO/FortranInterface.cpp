#include "CepGen/IO/FortranInterface.h"

#include "CepGen/StructureFunctions/StructureFunctionsBuilder.h"

#include "CepGen/StructureFunctions/SuriYennie.h"
#include "CepGen/StructureFunctions/SzczurekUleshchenko.h"
#include "CepGen/StructureFunctions/FioreBrasse.h"
#include "CepGen/StructureFunctions/ChristyBosted.h"
#include "CepGen/StructureFunctions/BlockDurandHa.h"
#include "CepGen/StructureFunctions/ALLM.h"

void
cepgen_structure_functions_( int& sfmode, double& q2, double& xbj, double& f2, double& fl )
{
  const CepGen::StructureFunctions sf = CepGen::StructureFunctionsBuilder::get( ( CepGen::StructureFunctionsType )sfmode, q2, xbj );
  f2 = sf.F2;
  fl = sf.FL;
}

