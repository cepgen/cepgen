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

#include "CepGen/Modules/CouplingFactory.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/Coupling.h"

namespace cepgen {
  class AlphaEMFixed final : public Coupling {
  public:
    explicit AlphaEMFixed(const ParametersList& params) : Coupling(params), value_(steer<double>("value")) {}

    static ParametersDescription description() {
      auto desc = Coupling::description();
      desc.setDescription("Constant alpha(EM)");
      desc.add<double>("value", constants::ALPHA_EM).setDescription("Constant value for alpha(EM)");
      return desc;
    }

    double operator()(double /* q */) const override { return value_; }

  private:
    double value_;
  };
}  // namespace cepgen

REGISTER_ALPHAEM_MODULE("fixed", AlphaEMFixed)
