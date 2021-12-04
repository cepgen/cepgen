/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

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

      static ParametersDescription parametersDescription();

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

    ParametersDescription BrashEtAl::parametersDescription() {
      auto desc = Parameterisation::parametersDescription();
      desc.setDescription("Brash et al.");
      return desc;
    }
  }  // namespace formfac
}  // namespace cepgen

REGISTER_FF_MODEL("Brash", BrashEtAl)
