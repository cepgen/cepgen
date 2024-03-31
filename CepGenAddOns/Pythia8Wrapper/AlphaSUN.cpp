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

#include <Pythia8/Pythia.h>
#include <Pythia8/StandardModel.h>

#include "CepGen/Modules/CouplingFactory.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/Coupling.h"

namespace cepgen {
  namespace pythia8 {
    class AlphaSUN final : public Coupling {
    public:
      explicit AlphaSUN(const ParametersList& params) : Coupling(params), alphas_(new Pythia8::AlphaSUN) {
        const auto nCHV = steer<int>("Ngauge"), alphaHVorder = nCHV > 1 ? steer<int>("alphaOrder") : 0;
        if (steer<bool>("setLambda")) {
          lambda_ = steer<double>("Lambda");
          alphas_->initLambda(nCHV, steer<int>("nFlav"), alphaHVorder, lambda_);
        } else {
          alphas_->initAlpha(
              nCHV, steer<int>("nFlav"), alphaHVorder, steer<double>("alphaFSR"), steer<double>("alphaFSRrefScale"));
          lambda_ = alphas_->Lambda();
        }
      }

      inline static ParametersDescription description() {
        auto desc = cepgen::Coupling::description();
        desc.setDescription("Pythia8 modelling of alpha(S) running in SU(N) model");
        desc.add<int>("Ngauge", 1);
        desc.add<int>("nFlav", 1);
        desc.add<int>("alphaOrder", 0);
        desc.add<bool>("setLambda", false);
        desc.add<double>("Lambda", 0.4);
        desc.add<double>("alphaFSR", 0.1);
        desc.add<double>("alphaFSRrefScale", 91.188);
        return desc;
      }

      inline double operator()(double q) const override { return alphas_->alpha(q * q); }

    private:
      const std::unique_ptr<Pythia8::AlphaSUN> alphas_;
      double lambda_{0.};
    };
  }  // namespace pythia8
}  // namespace cepgen
using Pythia8AlphaSUN = cepgen::pythia8::AlphaSUN;
REGISTER_ALPHAS_MODULE("pythia8UN", Pythia8AlphaSUN);
