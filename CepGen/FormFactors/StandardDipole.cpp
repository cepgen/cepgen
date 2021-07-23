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
        GE = pow(1. + q2 * inv_sq_scale_param_, -2.);
        GM = MU * GE;
      }

      static constexpr double inv_sq_scale_param_ = 1. / 0.71;  // for r_p = 0.81 fm
      //static constexpr double inv_sq_scale_param_ = 1. / 0.66; // for r_p = 0.84 fm
    };
  }  // namespace formfac
}  // namespace cepgen

REGISTER_FF_MODEL(gFFStandardDipoleHandler, StandardDipole)
