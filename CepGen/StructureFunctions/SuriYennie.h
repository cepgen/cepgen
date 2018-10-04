#ifndef CepGen_StructureFunctions_SuriYennie_h
#define CepGen_StructureFunctions_SuriYennie_h

#include "StructureFunctions.h"

namespace CepGen
{
  namespace SF
  {
    /// \f$F_{1/2/E/M}\f$ modelling by Suri and Yennie \cite Suri:1971yx
    class SuriYennie : public StructureFunctions
    {
      public:
        /// Collection of parameterisation-dependent couplings
        struct Parameterisation {
          // values extracted from experimental fits
          static Parameterisation standard();
          static Parameterisation alternative();

          double C1, C2, D1, rho2, Cp, Bp;
        };

        explicit SuriYennie( const SuriYennie::Parameterisation& param = SuriYennie::Parameterisation::standard() );
        SuriYennie& operator()( double xbj, double q2 ) override;

        double F1;
        double FE; ///< Electric proton form factor
        double FM; ///< Magnetic proton form factor
      private:
        Parameterisation params_;
    };
  }
}

#endif
