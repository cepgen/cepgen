/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2025  Laurent Forthomme
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

#include <boost/math/differentiation/finite_difference.hpp>

#include "CepGen/Modules/DerivatorFactory.h"
#include "CepGen/Utils/Derivator.h"
#include "CepGen/Utils/FunctionWrapper.h"

namespace cepgen::boost {
  /// Finite difference derivation algorithm
  /// \tparam N order of accuracy
  template <size_t N>
  class FiniteDifferenceDerivator final : public utils::Derivator {
  public:
    explicit FiniteDifferenceDerivator(const ParametersList& params) : utils::Derivator(params) {}

    static ParametersDescription description() {
      auto desc = utils::Derivator::description();
      desc.setDescription("Boost complex step derivation algorithm");
      return desc;
    }

    /// Evaluate the derivative of a function at a given value
    /// \param[in] function function to derive
    /// \param[in] x_coordinate coordinate
    /// \param[in] step_size (optional) step size ; if not provided, will use default algorithm value
    double derivate(const utils::FunctionWrapper& function, double x_coordinate, double /*step_size*/) const override {
      double uncertainty;
      return ::boost::math::differentiation::finite_difference_derivative<decltype(function), double, N>(
          function, x_coordinate, &uncertainty);
    }
  };
}  // namespace cepgen::boost
using FDDerivator1 = cepgen::boost::FiniteDifferenceDerivator<1>;
using FDDerivator2 = cepgen::boost::FiniteDifferenceDerivator<2>;
using FDDerivator4 = cepgen::boost::FiniteDifferenceDerivator<4>;
using FDDerivator6 = cepgen::boost::FiniteDifferenceDerivator<6>;
using FDDerivator8 = cepgen::boost::FiniteDifferenceDerivator<8>;
REGISTER_DERIVATOR("boost-finitediff1", FDDerivator1);
REGISTER_DERIVATOR("boost-finitediff2", FDDerivator2);
REGISTER_DERIVATOR("boost-finitediff4", FDDerivator4);
REGISTER_DERIVATOR("boost-finitediff6", FDDerivator6);
REGISTER_DERIVATOR("boost-finitediff8", FDDerivator8);
