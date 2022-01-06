/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2022  Laurent Forthomme
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
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/Utils.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

namespace cepgen {
  namespace strfun {
    /// \f$F_{1,2,E,M}\f$ modelling by Suri and Yennie \cite Suri:1971yx
    class SuriYennie : public Parameterisation {
    public:
      /// User-steered Suri-Yennie continuum structure functions calculator
      explicit SuriYennie(const ParametersList& params)
          : Parameterisation(params),
            C1_(steer<double>("C1")),
            C2_(steer<double>("C2")),
            D1_(steer<double>("D1")),
            rho2_(steer<double>("rho2")),
            Cp_(steer<double>("Cp")),
            Bp_(steer<double>("Bp")) {}

      static ParametersDescription description() {
        auto desc = Parameterisation::description();
        desc.setDescription("Suri-Yennie FE/FM");
        desc.add<double>("C1", 0.86926);
        desc.add<double>("C2", 2.23422);
        desc.add<double>("D1", 0.12549);
        desc.add<double>("rho2", 0.585);
        desc.add<double>("Cp", 0.96);
        desc.add<double>("Bp", 0.63);
        return desc;
      }

      SuriYennie& eval(double xbj, double q2) override {
        const double mx2 = utils::mX2(xbj, q2, mp2_), dm2 = mx2 - mp2_;  // [GeV^2]
        const double en = q2 + dm2;                                      // [GeV^2]
        const double nu = 0.5 * en / mp_, x_pr = q2 / (q2 + mx2), tau = 0.25 * q2 / mp2_;
        const double mq = rho2_ + q2;

        const double inv_q2 = 1. / q2;

        FM = inv_q2 *
             (C1_ * dm2 * pow(rho2_ / mq, 2) + C2_ * mp2_ * pow(1. - x_pr, 4) / (1. + x_pr * (x_pr * Cp_ - 2. * Bp_)));
        FE = (tau * FM + D1_ * dm2 * q2 * rho2_ / mp2_ * pow(dm2 / mq / en, 2)) / (1. + nu * nu * inv_q2);

        W1 = 0.5 * FM * q2 / mp_;
        W2 = 2. * mp_ * FE;
        F2 = 2. * nu * FE;
        return *this;
      }

    private:
      double C1_{0.}, C2_{0.};
      double D1_{0.};
      double rho2_{0.};
      double Cp_{0.}, Bp_{0.};
    };

    class SuriYennieAlt final : public SuriYennie {
    public:
      explicit SuriYennieAlt(const ParametersList& params) : SuriYennie(params) {}

      static ParametersDescription description() {
        auto desc = SuriYennie::description();
        desc.setDescription("Suri-Yennie alternative FE/FM parameterisation");
        desc.add<double>("C1", 0.6303);
        desc.add<double>("C2", 2.3049);
        desc.add<double>("D1", 0.04681);
        desc.add<double>("rho2", 1.05);
        desc.add<double>("Cp", 1.23);
        desc.add<double>("Bp", 0.61);
        return desc;
      }
    };
  }  // namespace strfun
}  // namespace cepgen

REGISTER_STRFUN(strfun::Type::SuriYennie, SuriYennie, strfun::SuriYennie)
REGISTER_STRFUN(strfun::Type::SuriYennieAlt, SuriYennieAlt, strfun::SuriYennieAlt)
