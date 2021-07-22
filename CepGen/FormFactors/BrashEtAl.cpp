#include <cmath>

#include "CepGen/Core/Exception.h"
#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"

namespace cepgen {
  namespace formfac {
    /// \cite Brash:2001qq
    class BrashEtAl final : public Parameterisation {
    public:
      using Parameterisation::Parameterisation;
      static std::string description() { return "Brash et al."; }

    private:
      static const float MAX_Q2;
      void compute(double q2) override {
        if (q2 > MAX_Q2)
          CG_WARNING("BrashEtAl") << "Q² = " << q2 << " > " << MAX_Q2 << " GeV² = max(Q²).\n\t"
                                  << "Brash et al. FF parameterisation not designed for high-Q² values.";
        const double q = sqrt(q2);
        GM = 1. / (1. + q * (0.116 + q * (2.874 + q * (0.241 + q * (1.006 + q * 0.345)))));

        const double r = std::min(1., 1. - 0.13 * (q2 - 0.04));
        if (r < 0.) {
          GM = GE = 0.;
          return;
        }
        GE = r * GM;
        GM *= MU;
      }
    };

    const float BrashEtAl::MAX_Q2 = 7.7;
  }  // namespace formfac
}  // namespace cepgen

REGISTER_FF_MODEL("Brash", BrashEtAl)
