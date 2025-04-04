/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2021-2025  Laurent Forthomme
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

#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h>

#include "CepGen/Modules/DerivatorFactory.h"
#include "CepGen/Utils/Derivator.h"
#include "CepGen/Utils/GSLFunctionsWrappers.h"
#include "CepGen/Utils/Message.h"

using namespace cepgen;
using namespace cepgen::utils;

class GSLDerivator final : public Derivator {
public:
  explicit GSLDerivator(const ParametersList& params) : Derivator(params), mode_(steerAs<int, Mode>("mode")) {}

  static ParametersDescription description() {
    auto desc = Derivator::description();
    desc.setDescription("GSL numerical differentiation algorithm");
    desc.addAs<int, Mode>("mode", Mode::central)
        .allow(0, "central")
        .allow(1, "forward")
        .allow(2, "backward")
        .setDescription("mode used for the adaptive difference algorithm");
    return desc;
  }

  /// Evaluate the derivative of a function at a given value
  /// \param[in] func function to derive
  /// \param[in] x coordinate
  /// \param[in] h (optional) step size ; if not provided, will use default algorithm value
  double derivate(const FunctionWrapper& func, double x, double h = -1.) const override {
    int res{GSL_SUCCESS};
    double val, val_unc;
    const double step_size = h > 0. ? h : h_;
    auto gfunc = GSLFunctionWrapper::build(func);
    switch (mode_) {
      case Mode::central:
        res = gsl_deriv_central(gfunc.get(), x, step_size, &val, &val_unc);
        break;
      case Mode::forward:
        res = gsl_deriv_forward(gfunc.get(), x, step_size, &val, &val_unc);
        break;
      case Mode::backward:
        res = gsl_deriv_backward(gfunc.get(), x, step_size, &val, &val_unc);
        break;
    }
    if (res != GSL_SUCCESS)
      CG_WARNING("GSLDerivator") << "Failed to evaluate the derivative. GSL error: " << gsl_strerror(res) << ".";
    return val;
  }

  enum struct Mode {
    central,  ///< adaptive central difference algorithm
    forward,  ///< adaptive forward difference algorithm
    backward  ///< adaptive backward difference algorithm
  };

private:
  const Mode mode_;
};
REGISTER_DERIVATOR("gsl", GSLDerivator);
