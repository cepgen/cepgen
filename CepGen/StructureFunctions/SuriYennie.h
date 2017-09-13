#ifndef CepGen_StructureFunctions_SuriYennie_h
#define CepGen_StructureFunctions_SuriYennie_h

#include "StructureFunctions.h"

namespace CepGen
{
  namespace SF
  {
    struct SuriYennieParameterisation {
      // values extracted from experimental fits
      static SuriYennieParameterisation standard();
      static SuriYennieParameterisation alternative();

      double C1, C2, D1, rho2, Cp, Bp;
    };

    StructureFunctions SuriYennie( double q2, double xbj, const SuriYennieParameterisation& param=SuriYennieParameterisation::standard() );
  }
}

#endif
