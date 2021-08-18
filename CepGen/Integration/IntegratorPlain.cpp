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

#include <gsl/gsl_monte_plain.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Integration/IntegratorGSL.h"
#include "CepGen/Modules/IntegratorFactory.h"

namespace cepgen {
  /// Plain integration algorithm randomly sampling points in the phase space
  class IntegratorPlain final : public IntegratorGSL {
  public:
    explicit IntegratorPlain(const ParametersList& params);
    static std::string description() { return "Plain (trial/error) integrator"; }

    void integrate(double& result, double& abserr) override;

  private:
    int ncvg_;
  };

  IntegratorPlain::IntegratorPlain(const ParametersList& params)
      : IntegratorGSL(params), ncvg_(params.get<int>("numFunctionCalls", 50000)) {}

  void IntegratorPlain::integrate(double& result, double& abserr) {
    //--- integration bounds
    std::vector<double> x_low(function_->dim, 0.), x_up(function_->dim, 1.);

    //--- launch integration
    std::unique_ptr<gsl_monte_plain_state, void (*)(gsl_monte_plain_state*)> pln_state(
        gsl_monte_plain_alloc(function_->dim), gsl_monte_plain_free);
    int res = gsl_monte_plain_integrate(
        function_.get(), &x_low[0], &x_up[0], function_->dim, ncvg_, gsl_rng_.get(), pln_state.get(), &result, &abserr);

    if (res != GSL_SUCCESS)
      throw CG_FATAL("Integrator:integrate") << "Error while performing the integration!\n\t"
                                             << "GSL error: " << gsl_strerror(res) << ".";

    result_ = result;
    err_result_ = abserr;
  }
}  // namespace cepgen

REGISTER_INTEGRATOR("plain", IntegratorPlain)
