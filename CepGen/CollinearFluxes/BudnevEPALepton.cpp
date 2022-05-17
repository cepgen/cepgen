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
#include "CepGen/Modules/CollinearFluxFactory.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/PDG.h"

namespace cepgen {
  namespace collflux {
    class BudnevEPALepton : public Parameterisation {
    public:
      explicit BudnevEPALepton(const ParametersList& params)
          : Parameterisation(params), ml2_(std::pow(PDG::get().mass(steer<int>("pdgId")), 2)) {
        CG_INFO("BudnevEPALepton") << "Budnev EPA for photon-from-lepton elastic limit (lepton: "
                                   << PDG::get().name(steer<int>("pdgId")) << ").\n\t "
                                   << "See V.M.Budnev, et al., Phys.Rep. 15C (1975) 181.";
      }

      static ParametersDescription description() {
        auto desc = Parameterisation::description();
        desc.setDescription("Budnev EPA for lepton");
        desc.add<int>("pdgId", 11).setDescription("lepton PDG id");
        return desc;
      }

      double operator()(double x, double = 0.) const override {
        if (x >= 1.)
          return 0.;
        double q2min = ml2_ * x * x / (1. - x);
        if (!t_range_.contains(q2min))
          return 0.;
        return std::max(0.,
                        0.5 * constants::ALPHA_EM * M_1_PI *
                            (2. * ml2_ * x * (-1. / q2min + 1. / t_range_.max()) +
                             (2. - 2. * x + x * x) / x * log(t_range_.max() / q2min)));
      }

    private:
      const double ml2_;
    };
  }  // namespace collflux
}  // namespace cepgen

REGISTER_COLLFLUX(BudnevEPALepton, collflux::BudnevEPALepton);
