/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022  Laurent Forthomme
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

#include "CepGen/CollinearFluxes/Parameterisation.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Modules/CollinearFluxFactory.h"
#include "CepGen/Physics/Beam.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/Limits.h"

namespace cepgen {
  namespace collflux {
    class BudnevEPA : public Parameterisation {
    public:
      explicit BudnevEPA(const ParametersList& params) : Parameterisation(params) {
        CG_INFO("BudnevEPA") << "Budnev EPA for photon-from-proton elastic limit.\n\t"
                             << "See V.M.Budnev, et al., Phys.Rep. 15C (1975) 181.";
      }

      static ParametersDescription description() {
        auto desc = Parameterisation::description();
        desc.setDescription("Budnev EPA for proton");
        return desc;
      }

      double operator()(double x, double = 0.) const override {
        if (x >= 1.)
          return 0.;
        double qmi = mp2_ * x * x / (1. - x);
        if (!t_range_.contains(qmi))
          return 0.;
        return std::max(
            0., constants::ALPHA_EM * M_1_PI * (phi_f(x, t_range_.max() / qz_) - phi_f(x, qmi / qz_)) * (1 - x) / x);
      }

    private:
      double phi_f(double x, double qq) const {
        const double qq1 = 1 + qq, y = x * x / (1 - x);
        double f = (1 + a_ * y) * (-log(qq1 / qq) + 1 / qq1 + 1 / (2 * qq1 * qq1) + 1 / (3 * qq1 * qq1 * qq1));
        f += (1 - b_) * y / (4 * qq * qq1 * qq1 * qq1);
        f += c_ * (1 + y / 4) *
             (log((qq1 - b_) / qq1) + b_ / qq1 + b_ * b_ / (2 * qq1 * qq1) + b_ * b_ * b_ / (3 * qq1 * qq1 * qq1));
        return f;
      }
      const double a_{7.16}, b_{-3.96}, c_{0.28};
      const double qz_{0.71};
    };
  }  // namespace collflux
}  // namespace cepgen

REGISTER_COLLFLUX(BudnevEPA, collflux::BudnevEPA);
