/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2022  Laurent Forthomme
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

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/GSLDerivator.h"
#include "CepGen/Utils/GSLFunctionsWrappers.h"

namespace cepgen {
  namespace utils {
    GSLDerivator::GSLDerivator(const ParametersList& params)
        : SteeredObject(params), mode_(steerAs<int, Mode>("mode")), h_(steer<double>("h")) {}

    ParametersDescription GSLDerivator::description() {
      auto desc = ParametersDescription();
      desc.add<int>("mode", (int)Mode::central);
      desc.add<double>("h", 1.e-2).setDescription("step size");
      return desc;
    }

    double GSLDerivator::eval(const Function1D& func, double x, double h) const {
      int res{GSL_SUCCESS};
      double val, val_unc;
      const double step_size = h > 0. ? h : h_;
      auto gfunc = utils::GSLFunctionWrapper::build(func);
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

  }  // namespace utils
}  // namespace cepgen
