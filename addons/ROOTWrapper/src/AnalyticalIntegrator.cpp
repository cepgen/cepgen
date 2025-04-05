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

#include <Math/Integrator.h>

#include "CepGen/Integration/AnalyticIntegrator.h"
#include "CepGen/Modules/AnalyticIntegratorFactory.h"
#include "CepGen/Utils/FunctionWrapper.h"
#include "CepGen/Utils/Message.h"

namespace cepgen::root {
  class AnalyticalIntegrator final : public AnalyticIntegrator {
  public:
    explicit AnalyticalIntegrator(const ParametersList& params)
        : AnalyticIntegrator(params),
          integrator_(steerAs<int, ROOT::Math::IntegrationOneDim::Type>("type"),
                      steer<double>("epsabs"),
                      steer<double>("epsrel"),
                      steer<int>("limit"),
                      steer<int>("rule")) {
      CG_DEBUG("root:AnalyticalIntegrator").log([this](auto& log) {
        log << "ROOT analytical integrator built with options:\n";
        integrator_.Options().Print(log.stream());
      });
    }

    static ParametersDescription description() {
      auto desc = AnalyticIntegrator::description();
      desc.setDescription("ROOT integration algorithms wrapper");
      desc.addAs<int>("type", ROOT::Math::IntegrationOneDim::Type::kDEFAULT).setDescription("type of integration");
      desc.add("epsabs", -1.).setDescription("desired absolute error limit");
      desc.add("epsrel", -1.).setDescription("desired relative error limit");
      desc.add("limit", 0).setDescription("maximum number of sub-intervals to build");
      desc.add("rule", 0).setDescription("Gauss-Kronrod integration rule (only for GSL kADAPTIVE type)");
      return desc;
    }

  private:
    double run(const utils::FunctionWrapper& integrand, void* params, const Limits& lim) const override {
      const auto func_local =
          utils::FunctionWrapper([&integrand, &params](double x) -> double { return integrand(x, params); });
      return integrator_.Integral(
          func_local, lim.hasMin() ? lim.min() : range_.min(), lim.hasMax() ? lim.max() : range_.max());
    }

    mutable ROOT::Math::IntegratorOneDim integrator_;
  };
}  // namespace cepgen::root
using cepgen::root::AnalyticalIntegrator;
REGISTER_ANALYTIC_INTEGRATOR("root", AnalyticalIntegrator);
