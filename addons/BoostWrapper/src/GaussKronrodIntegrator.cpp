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

#include "CepGen/Core/Exception.h"
#include "CepGen/Integration/BaseIntegrator.h"
#include "CepGen/Integration/Integrand.h"
#include "CepGen/Modules/BaseIntegratorFactory.h"
#include "CepGen/Utils/FunctionWrapper.h"
#include "CepGen/Utils/RandomGenerator.h"

namespace cepgen::boost {
  /// Boost Gauss-Kronrod integration algorithm
  template <size_t N>
  class GaussKronrodIntegrator final : public BaseIntegrator {
  public:
    explicit GaussKronrodIntegrator(const ParametersList& params)
        : BaseIntegrator(params), max_refinements_(steerAs<int, size_t>("limit")), tol_(steer<double>("tolerance")) {}

    static ParametersDescription description() {
      auto desc = BaseIntegrator::description();
      desc.setDescription("Boost Gauss-Kronrod integration algorithm");
      desc.add("limit", 1000).setDescription("maximum number of sub-intervals to build");
      desc.add("tolerance", std::numeric_limits<double>::infinity()).setDescription("maximal tolerance");
      return desc;
    }

  private:
    Value run(Integrand& integrand, const std::vector<Limits>& range = {}) const override {
      if (integrand.size() != 1)
        throw CG_FATAL("GaussKronrodIntegrator") << "This integration algorithm only runs on 1-dimensional integrands.";
      return Value{::boost::math::quadrature::gauss_kronrod<double, N>::integrate(
          [&integrand](double x) { return integrand.eval(std::vector{x}); },
          range.at(0).min(),
          range.at(0).max(),
          tol_,
          max_refinements_)};
    }
    const size_t max_refinements_;
    const double tol_;
  };
}  // namespace cepgen::boost
using BGKIntegrator15 = cepgen::boost::GaussKronrodIntegrator<15>;
using BGKIntegrator31 = cepgen::boost::GaussKronrodIntegrator<31>;
using BGKIntegrator41 = cepgen::boost::GaussKronrodIntegrator<41>;
using BGKIntegrator51 = cepgen::boost::GaussKronrodIntegrator<51>;
using BGKIntegrator61 = cepgen::boost::GaussKronrodIntegrator<61>;
REGISTER_BASE_INTEGRATOR("boost_gk15", BGKIntegrator15);
REGISTER_BASE_INTEGRATOR("boost_gk31", BGKIntegrator31);
REGISTER_BASE_INTEGRATOR("boost_gk41", BGKIntegrator41);
REGISTER_BASE_INTEGRATOR("boost_gk51", BGKIntegrator51);
REGISTER_BASE_INTEGRATOR("boost_gk61", BGKIntegrator61);
