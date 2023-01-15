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
#include "CepGen/Modules/FormFactorsFactory.h"

namespace cepgen {
  namespace formfac {
    /// \cite Brash:2001qq
    class BrashEtAl final : public Parameterisation {
    public:
      explicit BrashEtAl(const ParametersList& params)
          : Parameterisation(params),
            coeff_gm_(steer<std::vector<double> >("coeffGM")),
            coeff_r_(steer<std::vector<double> >("coeffR")),
            max_q2_(steer<double>("q2max")) {
        if (coeff_gm_.size() != 5)
          throw CG_FATAL("BrashEtAl") << "Invalid coefficients multiplicity for the G_M functional form!";
        if (coeff_r_.size() != 2)
          throw CG_FATAL("BrashEtAl") << "Invalid coefficients multiplicity for the G_E/G_M ratio functional form!";
      }

      static ParametersDescription description();

    private:
      FormFactors compute(double q2) override {
        if (q2 > max_q2_)
          CG_WARNING("BrashEtAl") << "Q² = " << q2 << " > " << max_q2_ << " GeV² = max(Q²).\n\t"
                                  << "Brash et al. FF parameterisation not designed for high-Q² values.";
        FormFactors out;
        const double r = std::min(1., 1. - coeff_r_.at(0) * (q2 - coeff_r_.at(1)));
        if (r < 0.)
          return out;
        const double q = sqrt(q2);
        out.GM =
            1. /
            (1. + q * (coeff_gm_.at(0) +
                       q * (coeff_gm_.at(1) + q * (coeff_gm_.at(2) + q * (coeff_gm_.at(3) + q * coeff_gm_.at(4))))));

        out.GE = r * out.GM;
        out.GM *= MU;
        return out;
      }
      const std::vector<double> coeff_gm_, coeff_r_;
      const double max_q2_;
    };

    ParametersDescription BrashEtAl::description() {
      auto desc = Parameterisation::description();
      desc.setDescription("Brash et al.");
      desc.add<std::vector<double> >("coeffGM", {0.116, 2.874, 0.241, 1.006, 0.345})
          .setDescription("coefficients for the G_M functional form");
      desc.add<std::vector<double> >("coeffR", {0.13, 0.04})
          .setDescription("coefficients for the G_E/G_M ratio functional form");
      desc.add<double>("q2max", 7.7).setDescription("maximal Q^2 supported (in GeV^2)");
      return desc;
    }
  }  // namespace formfac
}  // namespace cepgen

REGISTER_FORMFACTORS("Brash", BrashEtAl)
