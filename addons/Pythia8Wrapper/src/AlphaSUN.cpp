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
  class AlphaSUN final : public Coupling {
  public:
    explicit AlphaSUN(const ParametersList& params) : Coupling(params), alpha_s_(new Pythia8::AlphaSUN) {
      const auto nCHV = steer<int>("Ngauge"), alphaHVorder = nCHV > 1 ? steer<int>("alphaOrder") : 0;
      if (steer<bool>("setLambda")) {
        lambda_ = steer<double>("Lambda");
        alpha_s_->initLambda(nCHV, steer<int>("nFlav"), alphaHVorder, lambda_);
      } else {
        alpha_s_->initAlpha(
            nCHV, steer<int>("nFlav"), alphaHVorder, steer<double>("alphaFSR"), steer<double>("alphaFSRrefScale"));
        lambda_ = alpha_s_->Lambda();
      }
    }

    static ParametersDescription description() {
      auto desc = Coupling::description();
      desc.setDescription("Pythia8 modelling of alpha(S) running in SU(N) model");
      desc.add("Ngauge", 1);
      desc.add("nFlav", 1);
      desc.add("alphaOrder", 0);
      desc.add("setLambda", false);
      desc.add("Lambda", 0.4);
      desc.add("alphaFSR", 0.1);
      desc.add("alphaFSRrefScale", 91.188);
      return desc;
    }

    double operator()(double q) const override { return alpha_s_->alpha(q * q); }

  private:
    const std::unique_ptr<Pythia8::AlphaSUN> alpha_s_;
    double lambda_{0.};
  };
}  // namespace cepgen::pythia8
using Pythia8AlphaSUN = cepgen::pythia8::AlphaSUN;
REGISTER_ALPHAS_MODULE("pythia8UN", Pythia8AlphaSUN);
