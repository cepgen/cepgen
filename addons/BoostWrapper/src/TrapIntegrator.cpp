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

#include <boost/math/quadrature/trapezoidal.hpp>

#include "CepGen/Core/Exception.h"
#include "CepGen/Integration/BaseIntegrator.h"
#include "CepGen/Integration/Integrand.h"
#include "CepGen/Modules/BaseIntegratorFactory.h"
#include "CepGen/Utils/FunctionWrapper.h"
#include "CepGen/Utils/RandomGenerator.h"

namespace cepgen::boost {
  /// Trapezoidal integration algorithm
  class TrapIntegrator final : public BaseIntegrator {
  public:
    explicit TrapIntegrator(const ParametersList& params)
        : BaseIntegrator(params),
          max_refinements_(steerAs<int, size_t>("limit")),
          tolerance_(steer<double>("tolerance")) {}

    static ParametersDescription description() {
      auto desc = BaseIntegrator::description();
      desc.setDescription("Boost trapezoidal integration algorithm");
      desc.add("limit", 1000).setDescription("maximum number of sub-intervals to build");
      desc.add("tolerance", 1.e-6).setDescription("maximal tolerance");
      return desc;
    }

  private:
    Value run(Integrand& integrand, const std::vector<Limits>& range) const override {
      if (integrand.size() != 1)
        throw CG_FATAL("TrapIntegrator") << "This integration algorithm only runs on 1-dimensional integrands.";
      return Value{
          ::boost::math::quadrature::trapezoidal([&integrand](double x) { return integrand.eval(std::vector{x}); },
                                                 range.at(0).min(),
                                                 range.at(0).max(),
                                                 tolerance_,
                                                 max_refinements_)};
    }
    const size_t max_refinements_;
    const double tolerance_;
  };
}  // namespace cepgen::boost
using cepgen::boost::TrapIntegrator;
REGISTER_BASE_INTEGRATOR("boost", TrapIntegrator);
