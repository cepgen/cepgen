#ifndef CepGen_StructureFunctions_SuriYennie_h
#define CepGen_StructureFunctions_SuriYennie_h

#include "StructureFunctions.h"

namespace CepGen
{
  namespace SF
  {
    class SuriYennie : public StructureFunctions
    {
      public:
        struct Parameterisation {
          // values extracted from experimental fits
          static Parameterisation standard();
          static Parameterisation alternative();

          double C1, C2, D1, rho2, Cp, Bp;
        };

        SuriYennie( const SuriYennie::Parameterisation& param = SuriYennie::Parameterisation::standard() ) : params_( param ) {}
        SuriYennie operator()( double q2, double xbj ) const;

        double FE, FM;
      private:
        Parameterisation params_;
    };
  }
}

#endif
