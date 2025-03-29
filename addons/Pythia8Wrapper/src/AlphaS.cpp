/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024-2025  Laurent Forthomme
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

#include <Pythia8/Pythia.h>
#include <Pythia8/StandardModel.h>

#include "CepGen/Modules/CouplingFactory.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/Coupling.h"

namespace cepgen::pythia8 {
  class AlphaS final : public Coupling {
  public:
    explicit AlphaS(const ParametersList& params) : Coupling(params), alpha_s_(new Pythia8::AlphaStrong) {
      alpha_s_->init(
          steer<double>("alphaSvalue"), steer<int>("alphaSorder"), steer<int>("alphaSnfmax"), steer<bool>("useCMW"));
    }

    static ParametersDescription description() {
      auto desc = cepgen::Coupling::description();
      desc.setDescription("Pythia8 modelling of alpha(S) running");
      desc.add("alphaSvalue", 0.13);
      desc.add("alphaSorder", 1);
      desc.add("alphaSnfmax", 6);
      desc.add("useCMW", false);
      return desc;
    }

    double operator()(double q) const override { return alpha_s_->alphaS(q * q); }

  private:
    const std::unique_ptr<Pythia8::AlphaStrong> alpha_s_;
  };
}  // namespace cepgen::pythia8
using Pythia8AlphaS = cepgen::pythia8::AlphaS;
REGISTER_ALPHAS_MODULE("pythia8", Pythia8AlphaS);
