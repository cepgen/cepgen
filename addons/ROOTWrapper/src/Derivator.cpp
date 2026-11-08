/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2025  Laurent Forthomme
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

#include <TF1.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/DerivatorFactory.h"
#include "CepGen/Utils/Derivator.h"
#include "CepGen/Utils/FunctionWrapper.h"

namespace cepgen::root {
  class Derivator final : public utils::Derivator {
  public:
    explicit Derivator(const ParametersList& params) : utils::Derivator(params), order_(steer<int>("order")) {}

    static ParametersDescription description() {
      auto desc = utils::Derivator::description();
      desc.setDescription("ROOT derivation algorithm (Richardson's extrapolation method)");
      desc.add("order", 1).setDescription("order of the derivation");
      return desc;
    }

    /// Evaluate the derivative of a function at a given value
    /// \param[in] function Function to derive
    /// \param[in] x_coordinate Coordinate
    /// \param[in] step_size (Optional) step size; if not provided, will use default algorithm value
    Value derivate(const utils::FunctionWrapper& function, double x_coordinate, double step_size) const override {
      const auto root_function = TF1(
          "cepgen_functional",
          [&function](double vars[1], double* pars) { return function(vars[0], static_cast<void*>(pars)); },
          0.,
          1.,
          0);
      const auto epsilon = step_size < 0. ? h_ : step_size;
      double value;
      switch (order_) {
        case 1:
          value = root_function.Derivative(x_coordinate, nullptr, epsilon);
          break;
        case 2:
          value = root_function.Derivative2(x_coordinate, nullptr, epsilon);
          break;
        case 3:
          value = root_function.Derivative3(x_coordinate, nullptr, epsilon);
          break;
        default:
          throw CG_FATAL("root:Derivator") << "Invalid derivation order requested: " << order_ << ".";
      }
      return Value{value, root_function.DerivativeError()};
    }

  private:
    const int order_;
  };
}  // namespace cepgen::root
using ROOTDerivator = cepgen::root::Derivator;
REGISTER_DERIVATOR("root", ROOTDerivator);
