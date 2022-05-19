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
        : SteeredObject(params),
          range_(steer<Limits>("range")),
          mode_(steerAs<int, Mode>("mode")),
          fixed_type_(steerAs<int, FixedType>("fixedType")),
          limit_(steerAs<int, size_t>("limit")),
          epsabs_(steer<double>("epsabs")),
          epsrel_(steer<double>("epsrel")),
          func_params_(steer<ParametersList>("params")) {}

    ParametersDescription GSLIntegrator::description() {
      auto desc = ParametersDescription();
      desc.add<Limits>("range", Limits{0., 1.}).setDescription("integration range");
      desc.add<ParametersDescription>("params", ParametersDescription())
          .setDescription("parameters for the function to be integrated");
      desc.addAs<int, Mode>("mode", Mode::Fixed).setDescription("integrator algorithm to use");
      desc.addAs<int, FixedType>("fixedType", FixedType::Jacobi).setDescription("type of quadrature");
      desc.add<int>("limit", 1000).setDescription("maximum number of subintervals to build");
      desc.add<double>("epsabs", 0.).setDescription("desired absolute error limit");
      desc.add<double>("epsrel", 0.1).setDescription("desired relative error limit");
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
      int res = GSL_SUCCESS;
      if (mode_ == Mode::Fixed) {
        const gsl_integration_fixed_type* type{nullptr};
        switch (fixed_type_) {
          case FixedType::Legendre:
            type = gsl_integration_fixed_legendre;
            break;
          case FixedType::Chebyshev:
            type = gsl_integration_fixed_chebyshev;
            break;
          case FixedType::Gegenbauer:
            type = gsl_integration_fixed_gegenbauer;
            break;
          case FixedType::Jacobi:
            type = gsl_integration_fixed_jacobi;
            break;
          case FixedType::Laguerre:
            type = gsl_integration_fixed_laguerre;
            break;
          case FixedType::Hermite:
            type = gsl_integration_fixed_hermite;
            break;
          case FixedType::Exponential:
            type = gsl_integration_fixed_exponential;
            break;
          case FixedType::Rational:
            type = gsl_integration_fixed_rational;
            break;
          case FixedType::Chebyshev2:
            type = gsl_integration_fixed_chebyshev2;
            break;
          default:
            throw CG_FATAL("GSLIntegrator") << "Invalid fixed quadrature type: " << (int)fixed_type_ << ".";
        }
        std::unique_ptr<gsl_integration_fixed_workspace, void (*)(gsl_integration_fixed_workspace*)> workspace(
            gsl_integration_fixed_alloc(type, 50, xmin, xmax, 0., 0.), gsl_integration_fixed_free);
        res = gsl_integration_fixed(wrp, &result, workspace.get());
      } else if (mode_ == Mode::QNG) {
        size_t neval;
        double error;
        res = gsl_integration_qng(wrp, xmin, xmax, epsabs_, epsrel_, &result, &error, &neval);
      } else {
        double error;
        std::unique_ptr<gsl_integration_workspace, void (*)(gsl_integration_workspace*)> workspace(
            gsl_integration_workspace_alloc(limit_), gsl_integration_workspace_free);
        if (mode_ == Mode::QAG) {
          int key = GSL_INTEG_GAUSS41;
          res = gsl_integration_qag(wrp, xmin, xmax, epsabs_, epsrel_, limit_, key, workspace.get(), &result, &error);
        } else if (mode_ == Mode::QAGS)
          res = gsl_integration_qags(wrp, xmin, xmax, epsabs_, epsrel_, limit_, workspace.get(), &result, &error);
        else if (mode_ == Mode::QAWC)
          res = gsl_integration_qawc(
              const_cast<gsl_function*>(wrp), xmin, xmax, epsabs_, epsrel_, 0., limit_, workspace.get(), &result, &error);
      }
      if (res != GSL_SUCCESS)
        CG_WARNING("GSLIntegrator") << "Failed to evaluate the integral. GSL error: " << gsl_strerror(res) << ".";
      return result;
    }
  }  // namespace utils
}  // namespace cepgen
