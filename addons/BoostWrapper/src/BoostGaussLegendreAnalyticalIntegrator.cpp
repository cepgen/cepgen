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

#include <boost/math/quadrature/gauss.hpp>

#include "CepGen/Integration/AnalyticIntegrator.h"
#include "CepGen/Modules/AnalyticIntegratorFactory.h"
#include "CepGen/Utils/FunctionWrapper.h"

using namespace cepgen;

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

private:
  double run(const utils::FunctionWrapper& func, void* = nullptr, const Limits& lim = {}) const override {
    return boost::math::quadrature::gauss<double, N>::integrate(
        func, lim.hasMin() ? lim.min() : range_.min(), lim.hasMax() ? lim.max() : range_.max());
  }
};
using BGLIntegrator7 = BoostGaussLegendreAnalyticalIntegrator<7>;
using BGLIntegrator15 = BoostGaussLegendreAnalyticalIntegrator<15>;
using BGLIntegrator20 = BoostGaussLegendreAnalyticalIntegrator<20>;
using BGLIntegrator25 = BoostGaussLegendreAnalyticalIntegrator<25>;
using BGLIntegrator30 = BoostGaussLegendreAnalyticalIntegrator<30>;
REGISTER_ANALYTIC_INTEGRATOR("boost_gl7", BGLIntegrator7);
REGISTER_ANALYTIC_INTEGRATOR("boost_gl15", BGLIntegrator15);
REGISTER_ANALYTIC_INTEGRATOR("boost_gl20", BGLIntegrator20);
REGISTER_ANALYTIC_INTEGRATOR("boost_gl25", BGLIntegrator25);
REGISTER_ANALYTIC_INTEGRATOR("boost_gl30", BGLIntegrator30);
