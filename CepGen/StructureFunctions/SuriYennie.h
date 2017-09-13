#ifndef CepGen_StructureFunctions_SuriYennie_h
#define CepGen_StructureFunctions_SuriYennie_h

#include "StructureFunctions.h"

namespace CepGen
{
  namespace SF
  {
    struct SuriYennieParameterisation {
      double C1, C2, D1, rho2, Cp, Bp;
      // values extracted from experimental fits
      static SuriYennieParameterisation standard() {
        SuriYennieParameterisation p; p.C1 = 0.86926; p.C2 = 2.23422; p.D1 = 0.12549; p.rho2 = 0.585; p.Cp = 0.96; p.Bp = 0.63; return p;
      }
      static SuriYennieParameterisation alternative() {
        SuriYennieParameterisation p; p.C1 = 0.6303; p.C2 = 2.3049; p.D1 = 0.04681; p.rho2 = 1.05; p.Cp = 1.23; p.Bp = 0.61; return p;
      }
    };

    StructureFunctions SuriYennie( double q2, double xbj, const SuriYennieParameterisation& param=SuriYennieParameterisation::standard() );
  }
}

#endif
