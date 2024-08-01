/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2024  Laurent Forthomme
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

#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/Utils.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

namespace cepgen::strfun {
  /// \f$F_{1,2,E,M}\f$ modelling by Suri and Yennie \cite Suri:1971yx
  class SuriYennie : public Parameterisation {
  public:
    /// User-steered Suri-Yennie continuum structure functions calculator
    explicit SuriYennie(const ParametersList& params)
        : Parameterisation(params),
          inv_mp_(1. / mp_),
          c1_(steer<double>("C1")),
          c2_(steer<double>("C2")),
          d1_(steer<double>("D1")),
          rho2_(steer<double>("rho2")),
          cp_(steer<double>("Cp")),
          bp_(steer<double>("Bp")) {}

    static ParametersDescription description() {
      auto desc = Parameterisation::description();
      desc.setDescription("Suri-Yennie");
      desc.add<double>("C1", 0.86926);
      desc.add<double>("C2", 2.23422);
      desc.add<double>("D1", 0.12549);
      desc.add<double>("rho2", 0.585);
      desc.add<double>("Cp", 0.96);
      desc.add<double>("Bp", 0.63);
      return desc;
    }

    inline void eval() override {
      const double mx2 = utils::mX2(args_.xbj, args_.q2, mp2_), dm2 = mx2 - mp2_;  // [GeV^2]
      const double en = args_.q2 + dm2;                                            // [GeV^2]
      const double nu = 0.5 * en * inv_mp_, x_pr = args_.q2 / (args_.q2 + mx2),
                   tau = 0.25 * args_.q2 * inv_mp_ * inv_mp_;
      const double mq = rho2_ + args_.q2;

      const double inv_q2 = 1. / args_.q2;

      const double fm = inv_q2 * (c1_ * dm2 * std::pow(rho2_ / mq, 2) +
                                  c2_ * mp2_ * std::pow(1. - x_pr, 4) / (1. + x_pr * (x_pr * cp_ - 2. * bp_)));
      const double fe =
          (tau * fm + d1_ * dm2 * args_.q2 * rho2_ * std::pow(dm2 * inv_mp_ / mq / en, 2)) / (1. + nu * nu * inv_q2);

      setFE(fe);
      setFM(fm);
      setW1(0.5 * fm * args_.q2 * inv_mp_);
      setW2(2. * mp_ * fe);
      setF2(2. * nu * fe);
    }

  private:
    const double inv_mp_{0.};
    const double c1_{0.}, c2_{0.};
    const double d1_{0.};
    const double rho2_{0.};
    const double cp_{0.}, bp_{0.};
  };

  struct SuriYennieAlt final : public SuriYennie {
    using SuriYennie::SuriYennie;
    static ParametersDescription description() {
      auto desc = SuriYennie::description();
      desc.setDescription("Suri-Yennie (alternative)");
      desc.add<double>("C1", 0.6303);
      desc.add<double>("C2", 2.3049);
      desc.add<double>("D1", 0.04681);
      desc.add<double>("rho2", 1.05);
      desc.add<double>("Cp", 1.23);
      desc.add<double>("Bp", 0.61);
      return desc;
    }
  };
}  // namespace cepgen::strfun
using cepgen::strfun::SuriYennie;
using cepgen::strfun::SuriYennieAlt;
REGISTER_STRFUN("SuriYennie", 11, SuriYennie);
REGISTER_STRFUN("SuriYennieAlt", 14, SuriYennieAlt);
