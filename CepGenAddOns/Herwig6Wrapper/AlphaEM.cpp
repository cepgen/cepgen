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

#include "CepGen/Modules/CouplingFactory.h"
#include "CepGen/Physics/Coupling.h"
#include "CepGenAddOns/Herwig6Wrapper/Herwig6Interface.h"

namespace cepgen::herwig6 {
  class AlphaEM final : public cepgen::Coupling {
  public:
    explicit AlphaEM(const ParametersList& params) : cepgen::Coupling(params) {
      hwpram_.alphem = steer<double>("alphem");
    }

    inline static ParametersDescription description() {
      auto desc = cepgen::Coupling::description();
      desc.setDescription("Herwig6 modelling of alpha(EM) running");
      initialise();
      desc.add<double>("alphem", hwpram_.alphem).setDescription("alpha(EM) at beginning of evolution");
      return desc;
    }

    inline double operator()(double q) const override { return hwuaem(q * q); }
  };
}  // namespace cepgen::herwig6
using Herwig6AlphaEM = cepgen::herwig6::AlphaEM;
REGISTER_ALPHAEM_MODULE("herwig6", Herwig6AlphaEM);
