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

#include <cmath>

#include "CepGen/Modules/CouplingFactory.h"
#include "CepGen/Physics/Coupling.h"

namespace cepgen::apfelpp {
  class AlphaEM final : public cepgen::Coupling {
  public:
    explicit AlphaEM(const ParametersList& params)
        : cepgen::Coupling(params),
          alpha_em_(steer<double>("alphaQEDref"),
                    steer<double>("muQEDref"),
                    steer<std::vector<double> >("quarkThresholds"),
                    steer<std::vector<double> >("leptonThresholds"),
                    steer<int>("order")) {
      apfel::Banner();
    }

    static ParametersDescription description() {
      auto desc = cepgen::Coupling::description();
      desc.setDescription("APFEL++ alpha(EM) evolution algorithm");
      desc.add<double>("alphaQEDref", 1. / 128);
      desc.add<double>("muQEDref", 91.1876);
      desc.add<std::vector<double> >("quarkThresholds", {0., 0., 0., M_SQRT2, 4.5, 175.});
      desc.add<std::vector<double> >("leptonThresholds", {0., 0., 1.777});
      desc.add<int>("order", 0)
          .setDescription("QED evolution order")
          .allow(0, "leading order")
          .allow(1, "next-to-leading order");
      return desc;
    }

    inline double operator()(double q) const override { return alpha_em_.Evaluate(q); }

  private:
    const apfel::AlphaQED alpha_em_;
  };
}  // namespace cepgen::apfelpp
using AlphaEMAPFELpp = cepgen::apfelpp::AlphaEM;
REGISTER_ALPHAEM_MODULE("apfelpp", AlphaEMAPFELpp);
