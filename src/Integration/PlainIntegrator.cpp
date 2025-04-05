/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2025  Laurent Forthomme
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
#include "CepGen/Modules/BaseIntegratorFactory.h"
#include "CepGen/Utils/RandomGenerator.h"

using namespace cepgen;

/// Plain integration algorithm randomly sampling points in the phase space
class PlainIntegrator final : public GSLIntegrator {
public:
  explicit PlainIntegrator(const ParametersList& params)
      : GSLIntegrator(params), num_function_calls_(steer<int>("numFunctionCalls")) {}

  static ParametersDescription description() {
    auto desc = GSLIntegrator::description();
    desc.setDescription("Plain (trial/error) integrator");
    desc.add("numFunctionCalls", 50'000);
    return desc;
  }

  Value run(Integrand& integrand, const std::vector<Limits>& range) override {
    prepare(integrand, range);
    // launch integration
    const std::unique_ptr<gsl_monte_plain_state, decltype(&gsl_monte_plain_free)> plain_state(
        gsl_monte_plain_alloc(gsl_function_->dim), gsl_monte_plain_free);
    double result, absolute_error;
    if (const auto res = gsl_monte_plain_integrate(gsl_function_.get(),
                                                   &x_low_[0],
                                                   &x_high_[0],
                                                   gsl_function_->dim,
                                                   num_function_calls_,
                                                   random_generator_->engine<gsl_rng>(),
                                                   plain_state.get(),
                                                   &result,
                                                   &absolute_error);
        res != GSL_SUCCESS)
      throw CG_FATAL("Integrator:integrate") << "Error while performing the integration!\n\t"
                                             << "GSL error: " << gsl_strerror(res) << ".";
    return Value{result, absolute_error};
  }

private:
  const int num_function_calls_;
};
REGISTER_BASE_INTEGRATOR("plain", PlainIntegrator);
