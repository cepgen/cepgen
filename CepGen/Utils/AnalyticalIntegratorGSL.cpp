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
#include <gsl/gsl_version.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/SteeredObject.h"
#include "CepGen/Integration/AnalyticIntegrator.h"
#include "CepGen/Modules/AnalyticIntegratorFactory.h"
#include "CepGen/Utils/Functional.h"
#include "CepGen/Utils/GSLFunctionsWrappers.h"

namespace cepgen {
  namespace utils {
    class AnalyticalIntegratorGSL final : public AnalyticIntegrator {
    public:
      explicit AnalyticalIntegratorGSL(const ParametersList& params)
          : AnalyticIntegrator(params),
            mode_(steerAs<int, Mode>("mode")),
            fixed_type_(steerAs<int, FixedType>("fixedType")),
            nodes_(steer<int>("nodes")),
            alpha_(steer<double>("alpha")),
            beta_(steer<double>("beta")),
            limit_(steerAs<int, size_t>("limit")),
            epsabs_(steer<double>("epsabs")),
            epsrel_(steer<double>("epsrel")) {}

      static ParametersDescription description() {
        auto desc = AnalyticIntegrator::description();
        desc.setDescription("GSL 1D integration algorithms wrapper");
        desc.addAs<int, Mode>("mode", Mode::Fixed).setDescription("integrator algorithm to use");
        desc.addAs<int, FixedType>("fixedType", FixedType::Jacobi).setDescription("type of quadrature");
        desc.add<int>("nodes", 100).setDescription("number of quadrature nodes for the fixed type integration");
        desc.add<double>("alpha", 0.).setDescription("alpha parameter for the fixed type integration");
        desc.add<double>("beta", 0.).setDescription("alpha parameter for the fixed type integration");
        desc.add<int>("limit", 1000).setDescription("maximum number of subintervals to build");
        desc.add<double>("epsabs", 0.).setDescription("desired absolute error limit");
        desc.add<double>("epsrel", 0.1).setDescription("desired relative error limit");
        return desc;
      }

      double eval(const Function1D& func, void* obj = nullptr, const Limits& lim = {}) const override {
        if (obj)
          return eval(GSLFunctionWrapper::build(func, obj).get(), lim);
        return eval(GSLFunctionWrapper::build(func, func_params_).get(), lim);
      }

    private:
      enum struct Mode { Fixed = 0, QNG = 1, QAG = 2, QAGS = 3, QAWC = 4 };
      enum struct FixedType {
        Legendre = 0,
        Chebyshev = 1,
        Gegenbauer = 2,
        Jacobi = 3,
        Laguerre = 4,
        Hermite = 5,
        Exponential = 6,
        Rational = 7,
        Chebyshev2 = 8
      };
      double eval(const gsl_function*, const Limits&) const;
      const Mode mode_;
      const FixedType fixed_type_;
      const int nodes_;
      const double alpha_, beta_;
      const size_t limit_;
      const double epsabs_, epsrel_;
      static constexpr double INVALID = -999.999;
    };

    double AnalyticalIntegratorGSL::eval(const gsl_function* wrp, const Limits& lim) const {
      double result{0.};
#if defined(GSL_MAJOR_VERSION) && (GSL_MAJOR_VERSION > 2 || (GSL_MAJOR_VERSION == 2 && GSL_MINOR_VERSION >= 1))
      const double xmin = (lim.hasMin() ? lim.min() : range_.min());
      const double xmax = (lim.hasMax() ? lim.max() : range_.max());
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
            throw CG_FATAL("AnalyticalIntegratorGSL") << "Invalid fixed quadrature type: " << (int)fixed_type_ << ".";
        }
        std::unique_ptr<gsl_integration_fixed_workspace, void (*)(gsl_integration_fixed_workspace*)> workspace(
            gsl_integration_fixed_alloc(type, nodes_, xmin, xmax, alpha_, beta_), gsl_integration_fixed_free);
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
        CG_WARNING("AnalyticalIntegratorGSL")
            << "Failed to evaluate the integral. GSL error: " << gsl_strerror(res) << ".";
#else
      (void)lim;
      (void)wrp;
      CG_WARNING("AnalyticalIntegratorGSL") << "GSL version above 2.1 is required for integration.";
#endif
      return result;
    }
  }  // namespace utils
}  // namespace cepgen

REGISTER_ANALYTIC_INTEGRATOR("gsl", AnalyticalIntegratorGSL, utils::AnalyticalIntegratorGSL)
