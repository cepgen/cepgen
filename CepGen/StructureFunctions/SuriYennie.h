#ifndef CepGen_StructureFunctions_SuriYennie_h
#define CepGen_StructureFunctions_SuriYennie_h

#include "StructureFunctions.h"

namespace CepGen
{
  namespace SF
  {
    class SuriYennie
    {
      public:
        struct Parameterisation {
          // values extracted from experimental fits
          static Parameterisation standard();
          static Parameterisation alternative();

          double C1, C2, D1, rho2, Cp, Bp;
        };

        SuriYennie( const SuriYennie::Parameterisation& param = SuriYennie::Parameterisation::standard() ) : params_( param ) {}
        StructureFunctions operator()( double q2, double xbj ) const;

      private:
        Parameterisation params_;
    };
  }
}

#endif
