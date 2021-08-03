/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
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

#include "CepGen/Core/Exception.h"
#include "CepGen/Integration/Integrand.h"
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGen/Parameters.h"

namespace cepgen {
  /// Boost's Naive integration algorithm
  class IntegratorNaive final : public Integrator {
  public:
    explicit IntegratorNaive(const ParametersList&);
    static std::string description() { return "\"Naive\" Boost integrator"; }

    void integrate(double&, double&) override;

  private:
    std::function<double(const std::vector<double>&)> funct_;
    typedef boost::math::quadrature::naive_monte_carlo<double, decltype(funct_)> nmc_t;
    std::vector<std::pair<double, double> > bounds_;
    std::unique_ptr<nmc_t> mc_;
  };

  IntegratorNaive::IntegratorNaive(const ParametersList& params)
      : Integrator(params), funct_([=](const std::vector<double>& x) -> double { return integrand_->eval(x); }) {
    //--- a bit of printout for debugging
    CG_DEBUG("Integrator:build") << "Boost's Naive integrator built.";
  }

  void IntegratorNaive::integrate(double& result, double& abserr) {
    if (!initialised_) {
      bounds_ = std::vector<std::pair<double, double> >(integrand_->size(), {0., 1.});
      mc_.reset(new nmc_t(funct_, bounds_, 1.e-2, true, 1));
      initialised_ = true;
    }

    auto task = mc_->integrate();

    result_ = result = task.get();
    err_result_ = abserr = mc_->current_error_estimate();
  }
}  // namespace cepgen

REGISTER_INTEGRATOR("Naive", IntegratorNaive)
