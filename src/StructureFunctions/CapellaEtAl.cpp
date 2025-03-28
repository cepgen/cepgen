/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2025  Laurent Forthomme
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
#include "CepGen/StructureFunctions/Parameterisation.h"

namespace cepgen::strfun {
  /// Low-x structure functions, valid for the range \f$0<Q^2<5\f$ GeV^2
  /// \cite Capella:1994cr
  class CapellaEtAl final : public Parameterisation {
  public:
    explicit CapellaEtAl(const ParametersList& params)
        : Parameterisation(params),
          pA_(steer<double>("A")),
          pBu_(steer<double>("Bu")),
          pBd_(steer<double>("Bd")),
          alpha_r_(steer<double>("alphaR")),
          delta_0_(steer<double>("delta0")),
          coefficients_(steer<std::vector<double> >("coefficients")) {
      if (coefficients_.size() < 4)
        throw CG_FATAL("CapellaEtAl") << "Invalid multiplicity of coefficients given: " << coefficients_ << ".";
    }

    static ParametersDescription description() {
      auto desc = Parameterisation::description();
      desc.setDescription("Capella et al.");
      desc.add("A", 0.1502);
      desc.add("Bu", 1.2064);
      desc.add("Bd", 0.1798);
      desc.add("alphaR", 0.4150).setDescription("Reggeon intercept");
      desc.add("delta0", 0.08).setDescription("effective intercept at Q^2=0");
      desc.add("coefficients", std::vector{0.2631, 0.6452, 3.5489, 1.1170});
      return desc;
    }

    void eval() override {
      const auto nq2 = 1.5 * (1. + args_.q2 / (args_.q2 + coefficients_.at(2)));  // n(Q^2) function in the paper
      const auto dq2 = delta_0_ * (1 + (2 * args_.q2) / (args_.q2 + coefficients_.at(0)));  // big-Delta(Q^2) function
      const auto c1 = std::pow(args_.q2 / (args_.q2 + coefficients_.at(0)), 1. + dq2);
      const auto c2 = std::pow(args_.q2 / (args_.q2 + coefficients_.at(1)), alpha_r_);

      setF2(pA_ * std::pow(args_.xbj, -dq2) * std::pow(1 - args_.xbj, nq2 + 4.0) * c1 +
            std::pow(args_.xbj, 1.0 - alpha_r_) *
                (pBu_ * std::pow(1 - args_.xbj, nq2) + pBd_ * std::pow(1 - args_.xbj, nq2 + 1.0)) * c2);

      // Note: add the high-Q^2 QCD evolution case (eq. (8) in the paper)
    }

  private:
    const double pA_, pBu_, pBd_, alpha_r_, delta_0_;
    const std::vector<double> coefficients_;
  };
}  // namespace cepgen::strfun
using cepgen::strfun::CapellaEtAl;
REGISTER_STRFUN("Capella", 106, CapellaEtAl);
