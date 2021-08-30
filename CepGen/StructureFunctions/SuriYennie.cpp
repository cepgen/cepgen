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
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/StructureFunctions/SuriYennie.h"
#include "CepGen/Utils/Physics.h"

namespace cepgen {
  namespace strfun {
    SuriYennie::SuriYennie(const ParametersList &params) : Parameterisation(params) {
      const auto &model = params.get<std::string>("model", "standard");
      if (model == "standard") {
        sy_params_.C1 = 0.86926;
        sy_params_.C2 = 2.23422;
        sy_params_.D1 = 0.12549;
        sy_params_.rho2 = 0.585;
        sy_params_.Cp = 0.96;
        sy_params_.Bp = 0.63;
      } else if (model == "alternative") {
        sy_params_.C1 = 0.6303;
        sy_params_.C2 = 2.3049;
        sy_params_.D1 = 0.04681;
        sy_params_.rho2 = 1.05;
        sy_params_.Cp = 1.23;
        sy_params_.Bp = 0.61;
      } else {  // custom model
        sy_params_.C1 = params.get<double>("C1");
        sy_params_.C2 = params.get<double>("C2");
        sy_params_.D1 = params.get<double>("D1");
        sy_params_.rho2 = params.get<double>("rho2");
        sy_params_.Cp = params.get<double>("Cp");
        sy_params_.Bp = params.get<double>("Bp");
      }
    }

    SuriYennie &SuriYennie::eval(double xbj, double q2) {
      const double mx2 = utils::mX2(xbj, q2, mp2_), dm2 = mx2 - mp2_;  // [GeV^2]
      const double en = q2 + dm2;                                      // [GeV^2]
      const double nu = 0.5 * en / mp_, x_pr = q2 / (q2 + mx2), tau = 0.25 * q2 / mp2_;
      const double mq = sy_params_.rho2 + q2;

      const double inv_q2 = 1. / q2;

      FM = (sy_params_.C1 * dm2 * pow(sy_params_.rho2 / mq, 2) +
            (sy_params_.C2 * mp2_ * pow(1. - x_pr, 4)) / (1. + x_pr * (x_pr * sy_params_.Cp - 2. * sy_params_.Bp))) *
           inv_q2;
      FE = (tau * FM + sy_params_.D1 * dm2 * q2 * sy_params_.rho2 / mp2_ * pow(dm2 / mq / en, 2)) /
           (1. + nu * nu * inv_q2);

      W1 = 0.5 * FM * q2 / mp_;
      W2 = 2. * mp_ * FE;
      F2 = 2. * nu * FE;
      return *this;
    }
  }  // namespace strfun
}  // namespace cepgen

REGISTER_STRFUN(SuriYennie, strfun::SuriYennie)
