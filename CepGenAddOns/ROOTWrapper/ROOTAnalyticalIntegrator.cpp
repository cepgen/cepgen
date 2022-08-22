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

#include <Math/Integrator.h>

#include "CepGen/Integration/AnalyticIntegrator.h"
#include "CepGen/Modules/AnalyticIntegratorFactory.h"
#include "CepGen/Utils/FunctionsWrappers.h"

namespace cepgen {
  class ROOTAnalyticalIntegrator final : public AnalyticIntegrator {
  public:
    explicit ROOTAnalyticalIntegrator(const ParametersList& params)
        : AnalyticIntegrator(params),
          integr_(steerAs<int, ROOT::Math::IntegrationOneDim::Type>("type"),
                  steer<double>("epsabs"),
                  steer<double>("epsrel"),
                  steer<int>("limit")) {}

    static ParametersDescription description() {
      auto desc = AnalyticIntegrator::description();
      desc.setDescription("ROOT integration algorithms wrapper");
      desc.add<int>("type", -1).setDescription("type of integration");
      desc.add<int>("limit", 1000).setDescription("maximum number of subintervals to build");
      desc.add<double>("epsabs", 0.).setDescription("desired absolute error limit");
      desc.add<double>("epsrel", 0.1).setDescription("desired relative error limit");
      return desc;
    }

    double eval(const utils::Function1D& func, void* = nullptr, const Limits& lim = {}) const override {
      std::function<double(double)> func_local = func;
      const double xmin = (lim.hasMin() ? lim.min() : range_.min());
      const double xmax = (lim.hasMax() ? lim.max() : range_.max());
      return integr_.Integral(func_local, xmin, xmax);
    }

  private:
    mutable ROOT::Math::IntegratorOneDim integr_;
  };
}  // namespace cepgen

REGISTER_ANALYTIC_INTEGRATOR("root", ROOTAnalyticalIntegrator, ROOTAnalyticalIntegrator)
