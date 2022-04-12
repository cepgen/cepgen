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

#include <APFEL/APFEL.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/CouplingFactory.h"
#include "CepGen/Physics/Coupling.h"

namespace cepgen {
  class AlphaSAPFEL final : public Coupling {
  public:
    explicit AlphaSAPFEL(const ParametersList& params)
        : Coupling(params),
          order_(steer<int>("order", 2)),
          q0_(steer<double>("q0", 1.)),
          qmax_(steer<double>("qmax", 10000.)) {
      APFEL::SetPerturbativeOrder(order_);
      APFEL::InitializeAPFEL();
      APFEL::EvolveAPFEL(q0_, qmax_);
    }
    static ParametersDescription description() {
      auto desc = Coupling::description();
      desc.setDescription("APFEL alphaS evolution algorithm");
      desc.add<int>("order", 2).setDescription("perturbative evolution order");
      desc.add<double>("q0", 1.).setDescription("minimal Q reachable for evolution (in GeV)");
      desc.add<double>("qmax", 10000.).setDescription("maximal Q reachable (in GeV)");
      return desc;
    }

    double operator()(double q) const override {
      if (q < q0_ || q > qmax_)
        CG_WARNING("AlphaSAPFEL:get") << "q = " << q << " outside the evolution range"
                                      << " [" << q0_ << ":" << qmax_ << "].";
      return APFEL::AlphaQCD(q);
    }

  private:
    int order_;
    double q0_, qmax_;
  };
}  // namespace cepgen

REGISTER_ALPHAS_MODULE("apfel", AlphaSAPFEL)
