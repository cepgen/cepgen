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

#include <boost/math/quadrature/gauss.hpp>

#include "CepGen/Integration/AnalyticIntegrator.h"
#include "CepGen/Modules/AnalyticIntegratorFactory.h"
#include "CepGen/Utils/FunctionWrapper.h"

namespace cepgen {
  /// Gauss-Legendre integration algorithm
  template <size_t N>
  class BoostGaussLegendreAnalyticalIntegrator final : public AnalyticIntegrator {
  public:
    explicit BoostGaussLegendreAnalyticalIntegrator(const ParametersList& params) : AnalyticIntegrator(params) {}

    static ParametersDescription description() {
      auto desc = AnalyticIntegrator::description();
      desc.setDescription("Boost Gauss-Legendre integration algorithm");
      return desc;
    }

    double integrate(const utils::FunctionWrapper& func, void* = nullptr, const Limits& lim = {}) const override {
      return boost::math::quadrature::gauss<double, N>::integrate(
          func, lim.hasMin() ? lim.min() : range_.min(), lim.hasMax() ? lim.max() : range_.max());
    }
  };
}  // namespace cepgen
typedef cepgen::BoostGaussLegendreAnalyticalIntegrator<7> BGLIntegrator7;
typedef cepgen::BoostGaussLegendreAnalyticalIntegrator<15> BGLIntegrator15;
typedef cepgen::BoostGaussLegendreAnalyticalIntegrator<20> BGLIntegrator20;
typedef cepgen::BoostGaussLegendreAnalyticalIntegrator<25> BGLIntegrator25;
typedef cepgen::BoostGaussLegendreAnalyticalIntegrator<30> BGLIntegrator30;
REGISTER_ANALYTIC_INTEGRATOR("boost_gl7", BGLIntegrator7);
REGISTER_ANALYTIC_INTEGRATOR("boost_gl15", BGLIntegrator15);
REGISTER_ANALYTIC_INTEGRATOR("boost_gl20", BGLIntegrator20);
REGISTER_ANALYTIC_INTEGRATOR("boost_gl25", BGLIntegrator25);
REGISTER_ANALYTIC_INTEGRATOR("boost_gl30", BGLIntegrator30);
