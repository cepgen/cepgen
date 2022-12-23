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

#include <boost/math/quadrature/trapezoidal.hpp>

#include "CepGen/Integration/AnalyticIntegrator.h"
#include "CepGen/Modules/AnalyticIntegratorFactory.h"
#include "CepGen/Utils/FunctionsWrappers.h"

namespace cepgen {
  class BoostAnalyticalIntegrator final : public AnalyticIntegrator {
  public:
    explicit BoostAnalyticalIntegrator(const ParametersList& params)
        : AnalyticIntegrator(params),
          max_refinements_(steerAs<int, size_t>("limit")),
          tol_(steer<double>("tolerance")) {}

    static ParametersDescription description() {
      auto desc = AnalyticIntegrator::description();
      desc.setDescription("Boost trapezoidal integration algorithm");
      desc.add<int>("limit", 1000).setDescription("maximum number of subintervals to build");
      desc.add<double>("tolerance", 1.e-6).setDescription("maximal tolerance");
      return desc;
    }

    double integrate(const utils::Function1D& func, void* = nullptr, const Limits& lim = {}) const override {
      const double xmin = (lim.hasMin() ? lim.min() : range_.min());
      const double xmax = (lim.hasMax() ? lim.max() : range_.max());
      return boost::math::quadrature::trapezoidal(func, xmin, xmax, tol_, max_refinements_);
    }

  private:
    const size_t max_refinements_;
    const double tol_;
  };
}  // namespace cepgen

REGISTER_ANALYTIC_INTEGRATOR("boost", BoostAnalyticalIntegrator, BoostAnalyticalIntegrator)
