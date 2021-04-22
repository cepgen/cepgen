#ifndef CepGen_StructureFunctions_SuriYennie_h
#define CepGen_StructureFunctions_SuriYennie_h

#include "CepGen/StructureFunctions/Parameterisation.h"

namespace cepgen {
  namespace strfun {
    /// \f$F_{1,2,E,M}\f$ modelling by Suri and Yennie \cite Suri:1971yx
    class SuriYennie : public Parameterisation {
    public:
      /// Collection of parameterisation-dependent couplings
      struct Parameters {
        double C1, C2;
        double D1;
        double rho2;
        double Cp, Bp;
      };

      /// User-steered Suri-Yennie continuum structure functions calculator
      explicit SuriYennie(const ParametersList& params = ParametersList());
      static std::string description() { return "Suri-Yennie FE/FM"; }

      SuriYennie& eval(double xbj, double q2) override;

      double W1;  ///< Longitudinal form factor
      double W2;
      double FE;  ///< Electric proton form factor
      double FM;  ///< Magnetic proton form factor
    private:
      Parameters params_;
    };
  }  // namespace strfun
}  // namespace cepgen

#endif
