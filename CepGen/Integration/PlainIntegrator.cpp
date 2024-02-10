/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2023  Laurent Forthomme
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
#include "CepGen/Integration/GSLIntegrator.h"
#include "CepGen/Modules/IntegratorFactory.h"

namespace cepgen {
  /// Plain integration algorithm randomly sampling points in the phase space
  class PlainIntegrator final : public GSLIntegrator {
  public:
    explicit PlainIntegrator(const ParametersList& params)
        : GSLIntegrator(params), ncvg_(steer<int>("numFunctionCalls")) {}

    static ParametersDescription description() {
      auto desc = GSLIntegrator::description();
      desc.setDescription("Plain (trial/error) integrator");
      desc.add<int>("numFunctionCalls", 50'000);
      return desc;
    }

    Value integrate(Integrand& integrand) override {
      setIntegrand(integrand);

      //--- launch integration
      std::unique_ptr<gsl_monte_plain_state, decltype(&gsl_monte_plain_free)> pln_state(
          gsl_monte_plain_alloc(function_->dim), gsl_monte_plain_free);
      double result, abserr;
      if (int res = gsl_monte_plain_integrate(function_.get(),
                                              &xlow_[0],
                                              &xhigh_[0],
                                              function_->dim,
                                              ncvg_,
                                              rnd_gen_->engine<gsl_rng>(),
                                              pln_state.get(),
                                              &result,
                                              &abserr);
          res != GSL_SUCCESS)
        throw CG_FATAL("Integrator:integrate") << "Error while performing the integration!\n\t"
                                               << "GSL error: " << gsl_strerror(res) << ".";

      return Value{result, abserr};
    }

  private:
    const int ncvg_;
  };
}  // namespace cepgen
REGISTER_INTEGRATOR("plain", PlainIntegrator);
