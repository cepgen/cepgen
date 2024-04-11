/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024  Laurent Forthomme
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

#include <apfel/apfelxx.h>

#include "CepGen/Modules/CouplingFactory.h"
#include "CepGen/Physics/Coupling.h"

namespace cepgen {
  namespace apfelpp {
    class AlphaS final : public cepgen::Coupling {
    public:
      explicit AlphaS(const ParametersList& params)
          : cepgen::Coupling(params),
            alpha_s_(steer<double>("alphaSref"),
                     steer<double>("muQCDref"),
                     steer<std::vector<double> >("quarkThresholds"),
                     steer<int>("order")) {}

      static ParametersDescription description() {
        auto desc = cepgen::Coupling::description();
        desc.setDescription("APFEL++ alpha(S) evolution algorithm");
        desc.add<double>("alphaSref", 0.118);
        desc.add<double>("muQCDref", 91.1876);
        desc.add<std::vector<double> >("quarkThresholds", {0., 0., 0., M_SQRT2, 4.5, 175.});
        desc.add<int>("order", 2)
            .setDescription("QCD perturbative evolution order")
            .allow(0, "LO")
            .allow(1, "NLO")
            .allow(2, "NNLO")
            .allow(3, "NNNLO");
        return desc;
      }

      double operator()(double q) const override { return alpha_s_.Evaluate(q); }

    private:
      const apfel::AlphaQCD alpha_s_;
    };
  }  // namespace apfelpp
}  // namespace cepgen
using AlphaSAPFELpp = cepgen::apfelpp::AlphaS;
REGISTER_ALPHAS_MODULE("apfelpp", AlphaSAPFELpp);
