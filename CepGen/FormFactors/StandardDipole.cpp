#include <cmath>

#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"

namespace cepgen {
  namespace formfac {
    class StandardDipole final : public Parameterisation {
    public:
      using Parameterisation::Parameterisation;
      static std::string description() { return "Standard dipole"; }

    private:
      void compute(double q2) override {
        GE = pow(1. + q2 / 0.71, -2.);
        GM = MU * GE;
      }
    };
  }  // namespace formfac
}  // namespace cepgen

REGISTER_FF_MODEL(gFFStandardDipoleHandler, StandardDipole)
