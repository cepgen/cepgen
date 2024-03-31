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
    class AlphaEM final : public Coupling {
    public:
      explicit AlphaEM(const ParametersList& params)
          : Coupling(params), pythia_(new Pythia8::Pythia), alphaem_(new Pythia8::AlphaEM) {
        pythia_->settings.parm("StandardModel:alphaEM0", steer<double>("alphaEM0"));
        pythia_->settings.parm("StandardModel:alphaEMmZ", steer<double>("alphaEMmZ"));
        alphaem_->init(steer<int>("order"), &pythia_->settings);
      }

      inline static ParametersDescription description() {
        auto desc = cepgen::Coupling::description();
        desc.setDescription("Pythia8 modelling of alpha(EM) running");
        desc.add<int>("order", 1);
        desc.add<double>("alphaEM0", 0.00729735);
        desc.add<double>("alphaEMmZ", 0.00781751);
        return desc;
      }

      inline double operator()(double q) const override { return alphaem_->alphaEM(q * q); }

    private:
      const std::unique_ptr<Pythia8::Pythia> pythia_;
      const std::unique_ptr<Pythia8::AlphaEM> alphaem_;
    };
  }  // namespace pythia8
}  // namespace cepgen
using PythiaAlphaEM = cepgen::pythia8::AlphaEM;
REGISTER_ALPHAEM_MODULE("pythia8", PythiaAlphaEM);
