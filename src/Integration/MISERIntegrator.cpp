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

#include <gsl/gsl_monte_miser.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Integration/GSLIntegrator.h"
#include "CepGen/Modules/IntegratorFactory.h"

using namespace cepgen;

/// MISER integration algorithm developed by W.H. Press and G.R. Farrar, as documented in \cite Press:1989vk.
class MISERIntegrator final : public GSLIntegrator {
public:
  explicit MISERIntegrator(const ParametersList& params)
      : GSLIntegrator(params), num_function_calls_(steer<int>("numFunctionCalls")) {}

  static ParametersDescription description() {
    auto desc = GSLIntegrator::description();
    desc.setDescription("MISER adaptive importance sampling integrator");
    desc.add("numFunctionCalls", 50'000).setDescription("Number of function calls per phase space point evaluation");
    desc.add("estimateFraction", 0.1);
    desc.add("minCalls", 16 * 10);
    desc.add("minCallsPerBisection", 32 * 16 * 10);
    desc.add("alpha", 2.);
    desc.add("dither", 0.1);
    return desc;
  }

  Value integrate(Integrand& integrand) override {
    setIntegrand(integrand);
    const std::unique_ptr<gsl_monte_miser_state, decltype(&gsl_monte_miser_free)> miser_state(
        gsl_monte_miser_alloc(gsl_function_->dim), gsl_monte_miser_free);
    miser_state->verbose = verbosity_;
    gsl_monte_miser_params_get(miser_state.get(), &miser_params_);
    miser_params_.estimate_frac = steer<double>("estimateFraction");
    miser_params_.min_calls = steer<int>("minCalls");
    miser_params_.min_calls_per_bisection = steer<int>("minCallsPerBisection");
    miser_params_.alpha = steer<double>("alpha");
    miser_params_.dither = steer<double>("dither");
    gsl_monte_miser_params_set(miser_state.get(), &miser_params_);

    CG_DEBUG("Integrator:build") << "MISER parameters:\n\t"
                                 << "Number of calls: " << miser_params_.min_calls << ", "
                                 << "per bisection: " << miser_params_.min_calls_per_bisection << ",\n\t"
                                 << "Estimate fraction: " << miser_params_.estimate_frac << ",\n\t"
                                 << "α-value: " << miser_params_.alpha << ",\n\t"
                                 << "Dither: " << miser_params_.dither << ".";

    // launch the full integration
    double result, absolute_error;
    if (const auto res = gsl_monte_miser_integrate(gsl_function_.get(),
                                                   &x_low_[0],
                                                   &x_high_[0],
                                                   gsl_function_->dim,
                                                   num_function_calls_,
                                                   random_number_generator_->engine<gsl_rng>(),
                                                   miser_state.get(),
                                                   &result,
                                                   &absolute_error);
        res != GSL_SUCCESS)
      throw CG_FATAL("Integrator:integrate") << "Error while performing the integration!\n\t"
                                             << "GSL error: " << gsl_strerror(res) << ".";

    return Value{result, absolute_error};
  }

private:
  const int num_function_calls_;
  gsl_monte_miser_params miser_params_{};
};
REGISTER_INTEGRATOR("MISER", MISERIntegrator);
