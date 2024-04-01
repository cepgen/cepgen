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

namespace {
  extern "C" {
  double hwuaem_(double&);
  }
}  // namespace

namespace cepgen {
  namespace herwig6 {
    class AlphaEM final : public Coupling {
    public:
      explicit AlphaEM(const ParametersList& params) : Coupling(params) {}

      inline static ParametersDescription description() {
        auto desc = cepgen::Coupling::description();
        desc.setDescription("Herwig6 modelling of alpha(EM) running");
        return desc;
      }

      inline double operator()(double q) const override {
        double q2 = q * q;
        return hwuaem_(q2);
      }
    };
  }  // namespace herwig6
}  // namespace cepgen
using Herwig6AlphaEM = cepgen::herwig6::AlphaEM;
REGISTER_ALPHAEM_MODULE("herwig6", Herwig6AlphaEM);
