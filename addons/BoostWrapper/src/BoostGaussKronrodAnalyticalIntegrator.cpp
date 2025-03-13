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

#include <boost/math/quadrature/gauss_kronrod.hpp>

#include "CepGen/Integration/AnalyticIntegrator.h"
#include "CepGen/Modules/AnalyticIntegratorFactory.h"
#include "CepGen/Utils/FunctionWrapper.h"

using namespace cepgen;

/// Boost Gauss-Kronrod integration algorithm
template <size_t N>
class BoostGaussKronrodAnalyticalIntegrator final : public AnalyticIntegrator {
public:
  explicit BoostGaussKronrodAnalyticalIntegrator(const ParametersList& params)
      : AnalyticIntegrator(params), max_refinements_(steerAs<int, size_t>("limit")), tol_(steer<double>("tolerance")) {}

  static ParametersDescription description() {
    auto desc = AnalyticIntegrator::description();
    desc.setDescription("Boost Gauss-Kronrod integration algorithm");
    desc.add<int>("limit", 1000).setDescription("maximum number of subintervals to build");
    desc.add<double>("tolerance", std::numeric_limits<double>::infinity()).setDescription("maximal tolerance");
    return desc;
  }

private:
  double run(const utils::FunctionWrapper& func, void* = nullptr, const Limits& lim = {}) const override {
    return boost::math::quadrature::gauss_kronrod<double, N>::integrate(
        func, lim.hasMin() ? lim.min() : range_.min(), lim.hasMax() ? lim.max() : range_.max(), tol_, max_refinements_);
  }

  const size_t max_refinements_;
  const double tol_;
};
using BGKIntegrator15 = BoostGaussKronrodAnalyticalIntegrator<15>;
using BGKIntegrator31 = BoostGaussKronrodAnalyticalIntegrator<31>;
using BGKIntegrator41 = BoostGaussKronrodAnalyticalIntegrator<41>;
using BGKIntegrator51 = BoostGaussKronrodAnalyticalIntegrator<51>;
using BGKIntegrator61 = BoostGaussKronrodAnalyticalIntegrator<61>;
REGISTER_ANALYTIC_INTEGRATOR("boost_gk15", BGKIntegrator15);
REGISTER_ANALYTIC_INTEGRATOR("boost_gk31", BGKIntegrator31);
REGISTER_ANALYTIC_INTEGRATOR("boost_gk41", BGKIntegrator41);
REGISTER_ANALYTIC_INTEGRATOR("boost_gk51", BGKIntegrator51);
REGISTER_ANALYTIC_INTEGRATOR("boost_gk61", BGKIntegrator61);
