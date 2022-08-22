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

#include <boost/math/quadrature/gauss.hpp>

#include "CepGen/Integration/AnalyticIntegrator.h"
#include "CepGen/Modules/AnalyticIntegratorFactory.h"
#include "CepGen/Utils/FunctionsWrappers.h"

namespace cepgen {
  template <size_t N>
  class BoostGaussLegendreAnalyticalIntegrator final : public AnalyticIntegrator {
  public:
    explicit BoostGaussLegendreAnalyticalIntegrator(const ParametersList& params) : AnalyticIntegrator(params) {}

    static ParametersDescription description() {
      auto desc = AnalyticIntegrator::description();
      desc.setDescription("Boost Gauss-Legendre integration algorithm");
      return desc;
    }

    double eval(const utils::Function1D& func, void* = nullptr, const Limits& lim = {}) const override {
      const double xmin = (lim.hasMin() ? lim.min() : range_.min());
      const double xmax = (lim.hasMax() ? lim.max() : range_.max());
      return boost::math::quadrature::gauss<double, N>::integrate(func, xmin, xmax);
    }
  };
}  // namespace cepgen

REGISTER_ANALYTIC_INTEGRATOR("boost-gl7",
                             BoostGaussLegendreAnalyticalIntegrator7,
                             BoostGaussLegendreAnalyticalIntegrator<7>)
REGISTER_ANALYTIC_INTEGRATOR("boost-gl15",
                             BoostGaussLegendreAnalyticalIntegrator15,
                             BoostGaussLegendreAnalyticalIntegrator<15>)
REGISTER_ANALYTIC_INTEGRATOR("boost-gl20",
                             BoostGaussLegendreAnalyticalIntegrator20,
                             BoostGaussLegendreAnalyticalIntegrator<20>)
REGISTER_ANALYTIC_INTEGRATOR("boost-gl25",
                             BoostGaussLegendreAnalyticalIntegrator25,
                             BoostGaussLegendreAnalyticalIntegrator<25>)
REGISTER_ANALYTIC_INTEGRATOR("boost-gl30",
                             BoostGaussLegendreAnalyticalIntegrator30,
                             BoostGaussLegendreAnalyticalIntegrator<30>)
