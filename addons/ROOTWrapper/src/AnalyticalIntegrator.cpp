/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2024  Laurent Forthomme
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

#include <Math/Integrator.h>

#include "CepGen/Integration/AnalyticIntegrator.h"
#include "CepGen/Modules/AnalyticIntegratorFactory.h"
#include "CepGen/Utils/FunctionWrapper.h"
#include "CepGen/Utils/Message.h"

namespace cepgen::root {
  class AnalyticalIntegrator final : public cepgen::AnalyticIntegrator {
  public:
    explicit AnalyticalIntegrator(const ParametersList& params)
        : cepgen::AnalyticIntegrator(params),
          integr_(steerAs<int, ROOT::Math::IntegrationOneDim::Type>("type"),
                  steer<double>("epsabs"),
                  steer<double>("epsrel"),
                  steer<int>("limit"),
                  steer<int>("rule")) {
      CG_DEBUG("root:AnalyticalIntegrator").log([this](auto& log) {
        log << "ROOT analytical integrator built with options:\n";
        integr_.Options().Print(log.stream());
      });
    }

    static ParametersDescription description() {
      auto desc = cepgen::AnalyticIntegrator::description();
      desc.setDescription("ROOT integration algorithms wrapper");
      desc.addAs<int>("type", ROOT::Math::IntegrationOneDim::Type::kDEFAULT).setDescription("type of integration");
      desc.add<double>("epsabs", -1.).setDescription("desired absolute error limit");
      desc.add<double>("epsrel", -1.).setDescription("desired relative error limit");
      desc.add<int>("limit", 0).setDescription("maximum number of subintervals to build");
      desc.add<int>("rule", 0).setDescription("Gauss-Kronrod integration rule (only for GSL kADAPTIVE type)");
      return desc;
    }

    double integrate(const utils::FunctionWrapper& func, void* params, const Limits& lim) const override {
      const auto func_local = utils::FunctionWrapper([&func, &params](double x) { return func(x, params); });
      const double xmin = lim.hasMin() ? lim.min() : range_.min();
      const double xmax = lim.hasMax() ? lim.max() : range_.max();
      return integr_.Integral(func_local, xmin, xmax);
    }

  private:
    mutable ROOT::Math::IntegratorOneDim integr_;
  };
}  // namespace cepgen::root
using cepgen::root::AnalyticalIntegrator;
REGISTER_ANALYTIC_INTEGRATOR("root", AnalyticalIntegrator);
