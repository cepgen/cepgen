/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022  Laurent Forthomme
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
#include <gsl/gsl_integration.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/Functional.h"
#include "CepGen/Utils/GSLIntegrator.h"

namespace cepgen {
  namespace utils {
    GSLIntegrator::GSLIntegrator(const ParametersList& params)
        : SteeredObject(params), range_(steer<Limits>("range")), func_params_(steer<ParametersList>("params")) {}

    ParametersDescription GSLIntegrator::description() {
      auto desc = ParametersDescription();
      desc.add<Limits>("range", Limits{0., 1.}).setDescription("integration range");
      desc.add<ParametersDescription>("params", ParametersDescription())
          .setDescription("parameters for the function to be integrated");
      return desc;
    }

    double GSLIntegrator::eval(const Function1D& func, double xmin, double xmax) const {
      return eval(GSLFunctionWrapper::build(func, func_params_).get(), xmin, xmax);
    }

    double GSLIntegrator::eval(const Function1D& func, void* obj, double xmin, double xmax) const {
      return eval(GSLFunctionWrapper::build(func, obj).get(), xmin, xmax);
    }

    double GSLIntegrator::eval(const gsl_function* wrp, double xmin, double xmax) const {
      if (xmin == INVALID)
        xmin = range_.min();
      if (xmax == INVALID)
        xmax = range_.max();
      double result{0.};
      std::unique_ptr<gsl_integration_fixed_workspace, void (*)(gsl_integration_fixed_workspace*)> workspace(
          gsl_integration_fixed_alloc(gsl_integration_fixed_jacobi, 50, xmin, xmax, 0., 0.),
          gsl_integration_fixed_free);
      int res = gsl_integration_fixed(wrp, &result, workspace.get());
      if (res != GSL_SUCCESS)
        CG_WARNING("GSLIntegrator") << "Failed to evaluate the integral. GSL error: " << gsl_strerror(res) << ".";
      return result;
    }
  }  // namespace utils
}  // namespace cepgen
