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
#include "CepGen/Utils/RandomGenerator.h"

using namespace cepgen;

/// MISER integration algorithm developed by W.H. Press and G.R. Farrar, as documented in \cite Press:1989vk.
class MISERIntegrator final : public GSLIntegrator {
public:
  explicit MISERIntegrator(const ParametersList& params)
      : GSLIntegrator(params), num_function_calls_(steer<int>("numFunctionCalls")) {}

  static ParametersDescription description() {
    auto desc = GSLIntegrator::description();
    desc.setDescription("MISER adaptive importance sampling integrator");
    desc.add("numFunctionCalls", 50'000).setDescription("number of function calls per phase space point evaluation");
    desc.add("estimateFraction", 0.1)
        .setDescription(
            "fraction of the currently available number of function calls allocated to estimating the variance at each "
            "recursive step");
    desc.add("minCalls", 16 * 10)
        .setDescription("minimum number of function calls required for each estimate of the variance");
    desc.add("minCallsPerBisection", 32 * 16 * 10)
        .setDescription("minimum number of function calls required to proceed with a bisection step");
    desc.add("alpha", 2.)
        .setDescription(
            "how the estimated variances for the two sub-regions of a bisection are combined when allocating points");
    desc.add("dither", 0.1)
        .setDescription(
            "size of the random fractional variation into each bisection, which can be used to break the symmetry of "
            "integrands which are concentrated near the exact center of the hypercubic integration region");
    return desc;
  }

  Value run(Integrand& integrand, const std::vector<Limits>& range) override {
    prepare(integrand, range);
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
                                                   x_low_.data(),
                                                   x_high_.data(),
                                                   gsl_function_->dim,
                                                   num_function_calls_,
                                                   random_generator_->engine<gsl_rng>(),
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
