#ifndef CepGen_StructureFunctions_SuriYennie_h
#define CepGen_StructureFunctions_SuriYennie_h

#include "CepGen/StructureFunctions/StructureFunctions.h"

namespace cepgen
{
  namespace sf
  {
    /// \f$F_{1/2/E/M}\f$ modelling by Suri and Yennie \cite Suri:1971yx
    class SuriYennie : public Parameterisation
    {
      public:
        /// Collection of parameterisation-dependent couplings
        struct Parameters {
          // values extracted from experimental fits
          static Parameters standard();
          static Parameters alternative();

          double C1, C2, D1, rho2, Cp, Bp;
        };

        explicit SuriYennie( const Parameters& param = Parameters::standard() );
        SuriYennie& operator()( double xbj, double q2 ) override;

        double F1;
        double FE; ///< Electric proton form factor
        double FM; ///< Magnetic proton form factor
      private:
        Parameters params_;
    };
  }
}

#endif
