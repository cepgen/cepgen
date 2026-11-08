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

#include "CepGen/Core/Exception.h"
#include "CepGen/Integration/Integrand.h"
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGen/Utils/FunctionWrapper.h"
#include "CepGen/Utils/RandomGenerator.h"

namespace cepgen::boost {
  /// Gauss-Legendre integration algorithm
  template <size_t N>
  class GaussLegendreIntegrator final : public Integrator {
  public:
    explicit GaussLegendreIntegrator(const ParametersList& params) : Integrator(params) {}

    static ParametersDescription description() {
      auto desc = Integrator::description();
      desc.setDescription("Boost Gauss-Legendre integration algorithm");
      return desc;
    }

  private:
    bool oneDimensional() const override { return true; }
    Value run(Integrand& integrand, const std::vector<Limits>& range = {}) override {
      if (integrand.size() != 1)
        throw CG_ERROR("GaussLegendreIntegrator")
            << "This integration algorithm only runs on 1-dimensional integrands.";
      return Value{::boost::math::quadrature::gauss<double, N>::integrate(
          [&integrand](double x) { return integrand.eval(std::vector{x}); }, range.at(0).min(), range.at(0).max())};
    }
  };
}  // namespace cepgen::boost
using BGLIntegrator7 = cepgen::boost::GaussLegendreIntegrator<7>;
using BGLIntegrator15 = cepgen::boost::GaussLegendreIntegrator<15>;
using BGLIntegrator20 = cepgen::boost::GaussLegendreIntegrator<20>;
using BGLIntegrator25 = cepgen::boost::GaussLegendreIntegrator<25>;
using BGLIntegrator30 = cepgen::boost::GaussLegendreIntegrator<30>;
REGISTER_INTEGRATOR("boost_gl7", BGLIntegrator7);
REGISTER_INTEGRATOR("boost_gl15", BGLIntegrator15);
REGISTER_INTEGRATOR("boost_gl20", BGLIntegrator20);
REGISTER_INTEGRATOR("boost_gl25", BGLIntegrator25);
REGISTER_INTEGRATOR("boost_gl30", BGLIntegrator30);
