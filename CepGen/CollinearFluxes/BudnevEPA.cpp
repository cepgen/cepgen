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
      explicit BudnevEPA(const ParametersList& params) : Parameterisation(params), q2max_(steer<double>("q2max")) {
        CG_INFO("BudnevEPA") << "Budnev EPA for photon-from-proton elastic limit.\n\t"
                             << "See V.M.Budnev, et al., Phys.Rep. 15C (1975) 181.";
      }

      static ParametersDescription description() {
        auto desc = Parameterisation::description();
        desc.setDescription("Budnev EPA for proton");
        desc.add<double>("q2max", 25.);
        return desc;
      }

      double operator()(double x, double = 0.) const override {
        //integer nb_proton(2), nb_neutron(2)
        //common/to_heavyion_pdg/ nb_proton, nb_neutron
        //double precision mass_ion(2)
        //common/to_heavyion_mass/mass_ion

        const double qz = 0.71;
        //if (nb_proton(beamid) != 1 || nb_neutron(beamid) != 0) {
        //  xin = mass_ion(beamid);
        //  alpha = alpha * nb_proton(beamid)
        //}

        if (x >= 1.)
          return 0.;
        double qmi = mp2_ * x * x / (1. - x);
        if (qmi >= q2max_)
          return 0.;
        return std::max(0., constants::ALPHA_EM * M_1_PI * (phi_f(x, q2max_ / qz) - phi_f(x, qmi / qz)) * (1 - x) / x);
      }

    private:
      double phi_f(double x, double qq) const {
        const double a = 7.16, b = -3.96, c = 0.28;
        double qq1 = 1 + qq, y = x * x / (1 - x);
        double f = (1 + a * y) * (-log(qq1 / qq) + 1 / qq1 + 1 / (2 * qq1 * qq1) + 1 / (3 * qq1 * qq1 * qq1));
        f += (1 - b) * y / (4 * qq * qq1 * qq1 * qq1);
        f += c * (1 + y / 4) *
             (log((qq1 - b) / qq1) + b / qq1 + b * b / (2 * qq1 * qq1) + b * b * b / (3 * qq1 * qq1 * qq1));
        return f;
      }

      const double q2max_;
    };
  }  // namespace collflux
}  // namespace cepgen

REGISTER_COLLFLUX(BudnevEPA, collflux::BudnevEPA);
