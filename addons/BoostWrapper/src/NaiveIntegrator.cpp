/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2020-2025  Laurent Forthomme
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

#include "CepGen/Integration/BaseIntegrator.h"
#include "CepGen/Integration/Integrand.h"
#include "CepGen/Modules/BaseIntegratorFactory.h"

namespace cepgen::boost {
  /// Boost's Naive integration algorithm
  class NaiveIntegrator final : public BaseIntegrator {
  public:
    explicit NaiveIntegrator(const ParametersList& params) : BaseIntegrator(params) {}

    static ParametersDescription description() {
      auto desc = BaseIntegrator::description();
      desc.setDescription("'Naive' Boost integrator");
      return desc;
    }

    Value run(Integrand& integrand, const std::vector<Limits>& range) const override {
      std::vector<std::pair<double, double> > bounds;
      std::transform(range.begin(), range.end(), std::back_inserter(bounds), [](const auto& lim) { return lim.raw(); });
      const auto funct = [&integrand](const std::vector<double>& coord) -> double { return integrand.eval(coord); };
      ::boost::math::quadrature::naive_monte_carlo<double, decltype(funct)> mc(funct, bounds, 1.e-2, true, 1);
      auto task = mc.integrate();
      return Value{task.get(), mc.current_error_estimate()};
    }
  };
}  // namespace cepgen::boost
using cepgen::boost::NaiveIntegrator;
REGISTER_BASE_INTEGRATOR("Naive", NaiveIntegrator);
