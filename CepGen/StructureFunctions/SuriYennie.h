#ifndef CepGen_StructureFunctions_SuriYennie_h
#define CepGen_StructureFunctions_SuriYennie_h

#include "CepGen/StructureFunctions/Parameterisation.h"

namespace cepgen
{
  namespace strfun
  {
    /// \f$F_{1,2,E,M}\f$ modelling by Suri and Yennie \cite Suri:1971yx
    class SuriYennie : public Parameterisation
    {
      public:
        /// Collection of parameterisation-dependent couplings
        struct Parameters {
          static Parameters standard(); ///< Standard parameterisation extracted from experimental fits
          static Parameters alternative(); ///< Alternative parameterisation extracted from experimental fits
          double C1, C2;
          double D1;
          double rho2;
          double Cp, Bp;
        };

        /// User-steered Suri-Yennie continuum structure functions calculator
        explicit SuriYennie( const ParametersList& params = ParametersList() );
        SuriYennie& operator()( double xbj, double q2 ) override;

        double F1; ///< Longitudinal form factor
        double FE; ///< Electric proton form factor
        double FM; ///< Magnetic proton form factor
      private:
        Parameters params_;
    };
  }
}

#endif
