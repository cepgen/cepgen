/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2020-2023  Laurent Forthomme
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

#include <boost/math/quadrature/naive_monte_carlo.hpp>

#include "CepGen/Integration/Integrand.h"
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Modules/IntegratorFactory.h"

namespace cepgen {
  /// Boost's Naive integration algorithm
  class IntegratorNaive final : public Integrator {
  public:
    explicit IntegratorNaive(const ParametersList& params) : Integrator(params) {}

    static ParametersDescription description() {
      auto desc = Integrator::description();
      desc.setDescription("\"Naive\" Boost integrator");
      return desc;
    }

    void setLimits(const std::vector<Limits>& lims) override {
      Integrator::setLimits(lims);
      bounds_.clear();
      std::transform(
          limits_.begin(), limits_.end(), std::back_inserter(bounds_), [](const auto& lim) { return lim.raw(); });
    }
    void integrate(Integrand& integrand, double& result, double& abserr) override {
      checkLimits(integrand);

      auto funct = [&](const std::vector<double>& coord) -> double { return integrand.eval(coord); };
      boost::math::quadrature::naive_monte_carlo<double, decltype(funct)> mc(funct, bounds_, 1.e-2, true, 1);
      auto task = mc.integrate();

      result_ = result = task.get();
      err_result_ = abserr = mc.current_error_estimate();
    }

  private:
    std::vector<std::pair<double, double> > bounds_;
  };
}  // namespace cepgen

REGISTER_INTEGRATOR("Naive", IntegratorNaive);
