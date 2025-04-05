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

#include <gsl/gsl_version.h>

#if defined(GSL_MAJOR_VERSION) && (GSL_MAJOR_VERSION > 2 || (GSL_MAJOR_VERSION == 2 && GSL_MINOR_VERSION >= 4))

#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Integration/BaseIntegrator.h"
#include "CepGen/Integration/Integrand.h"
#include "CepGen/Modules/BaseIntegratorFactory.h"
#include "CepGen/Utils/GSLFunctionsWrappers.h"
#include "CepGen/Utils/RandomGenerator.h"

using namespace cepgen;
using namespace cepgen::utils;

class AnalyticalIntegratorGSL final : public BaseIntegrator {
public:
  explicit AnalyticalIntegratorGSL(const ParametersList& params)
      : BaseIntegrator(params),
        mode_(steerAs<int, Mode>("mode")),
        fixed_type_(steerAs<int, FixedType>("fixedType")),
        nodes_(steer<int>("nodes")),
        alpha_(steer<double>("alpha")),
        beta_(steer<double>("beta")),
        limit_(steerAs<int, size_t>("limit")),
        epsabs_(steer<double>("epsabs")),
        epsrel_(steer<double>("epsrel")) {}

  static ParametersDescription description() {
    auto desc = BaseIntegrator::description();
    desc.setDescription("GSL 1D integration algorithms wrapper");
    desc.addAs<int>("mode", Mode::Fixed).setDescription("integrator algorithm to use");
    desc.addAs<int>("fixedType", FixedType::Jacobi).setDescription("type of quadrature");
    desc.add("nodes", 100).setDescription("number of quadrature nodes for the fixed type integration");
    desc.add("alpha", 0.).setDescription("alpha parameter for the fixed type integration");
    desc.add("beta", 0.).setDescription("alpha parameter for the fixed type integration");
    desc.add("limit", 1000).setDescription("maximum number of sub-intervals to build");
    desc.add("epsabs", 0.).setDescription("desired absolute error limit");
    desc.add("epsrel", 0.1).setDescription("desired relative error limit");
    return desc;
  }

private:
  Value run(Integrand& integrand, const std::vector<Limits>& range) const override {
    return eval(
        GSLFunctionWrapper::build(FunctionWrapper{[&integrand](double x) { return integrand.eval(std::vector{x}); }},
                                  integrand_parameters_)
            .get(),
        range.at(0));
  }

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
  Value eval(const gsl_function* wrp, const Limits& range) const {
#if defined(GSL_MAJOR_VERSION) && (GSL_MAJOR_VERSION > 2 || (GSL_MAJOR_VERSION == 2 && GSL_MINOR_VERSION >= 1))
    double result{0.}, error{0.};
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
          throw CG_FATAL("AnalyticalIntegratorGSL")
              << "Invalid fixed quadrature type: " << static_cast<int>(fixed_type_) << ".";
      }
      std::unique_ptr<gsl_integration_fixed_workspace, void (*)(gsl_integration_fixed_workspace*)> workspace(
          gsl_integration_fixed_alloc(type, nodes_, range.min(), range.max(), alpha_, beta_),
          gsl_integration_fixed_free);
      res = gsl_integration_fixed(wrp, &result, workspace.get());
    } else if (mode_ == Mode::QNG) {
      size_t neval;
      res = gsl_integration_qng(wrp, range.min(), range.max(), epsabs_, epsrel_, &result, &error, &neval);
    } else {
      std::unique_ptr<gsl_integration_workspace, void (*)(gsl_integration_workspace*)> workspace(
          gsl_integration_workspace_alloc(limit_), gsl_integration_workspace_free);
      if (mode_ == Mode::QAG) {
        int key = GSL_INTEG_GAUSS41;
        res = gsl_integration_qag(
            wrp, range.min(), range.max(), epsabs_, epsrel_, limit_, key, workspace.get(), &result, &error);
      } else if (mode_ == Mode::QAGS)
        res = gsl_integration_qags(
            wrp, range.min(), range.max(), epsabs_, epsrel_, limit_, workspace.get(), &result, &error);
      else if (mode_ == Mode::QAWC)
        res = gsl_integration_qawc(const_cast<gsl_function*>(wrp),
                                   range.min(),
                                   range.max(),
                                   epsabs_,
                                   epsrel_,
                                   0.,
                                   limit_,
                                   workspace.get(),
                                   &result,
                                   &error);
    }
    if (res != GSL_SUCCESS)
      CG_WARNING("AnalyticalIntegratorGSL")
          << "Failed to evaluate the integral. GSL error: " << gsl_strerror(res) << ".";
    return Value{result, error};
#else
    (void)range;
    (void)wrp;
    CG_WARNING("AnalyticalIntegratorGSL") << "GSL version above 2.1 is required for integration.";
    return {};
#endif
  }
  const Mode mode_;
  const FixedType fixed_type_;
  const int nodes_;
  const double alpha_, beta_;
  const size_t limit_;
  const double epsabs_, epsrel_;
};

REGISTER_BASE_INTEGRATOR("gsl", AnalyticalIntegratorGSL);
#endif
