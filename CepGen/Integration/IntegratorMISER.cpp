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

#include <gsl/gsl_monte_miser.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Integration/IntegratorGSL.h"
#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGen/Parameters.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  /// MISER integration algorithm developed by W.H. Press and G.R. Farrar, as documented in \cite Press:1989vk.
  class IntegratorMISER final : public IntegratorGSL {
  public:
    explicit IntegratorMISER(const ParametersList&);

    static ParametersDescription description();

    void integrate(Integrand&, double&, double&) override;

  private:
    int ncvg_;
    gsl_monte_miser_params miser_params_;
    /// A trivial deleter for the MISER integrator
    struct gsl_monte_miser_deleter {
      inline void operator()(gsl_monte_miser_state* state) { gsl_monte_miser_free(state); }
    };
    std::unique_ptr<gsl_monte_miser_state, gsl_monte_miser_deleter> miser_state_;
  };

  IntegratorMISER::IntegratorMISER(const ParametersList& params)
      : IntegratorGSL(params), ncvg_(steer<int>("numFunctionCalls")) {}

  void IntegratorMISER::integrate(Integrand& integrand, double& result, double& abserr) {
    setIntegrand(integrand);
    miser_state_.reset(gsl_monte_miser_alloc(function_->dim));
    miser_state_->verbose = verbosity_;
    gsl_monte_miser_params_get(miser_state_.get(), &miser_params_);
    miser_params_.estimate_frac = steer<double>("estimateFraction");
    miser_params_.min_calls = steer<int>("minCalls");
    miser_params_.min_calls_per_bisection = steer<int>("minCallsPerBisection");
    miser_params_.alpha = steer<double>("alpha");
    miser_params_.dither = steer<double>("dither");
    gsl_monte_miser_params_set(miser_state_.get(), &miser_params_);

    CG_DEBUG("Integrator:build") << "MISER parameters:\n\t"
                                 << "Number of calls: " << miser_params_.min_calls << ", "
                                 << "per bisection: " << miser_params_.min_calls_per_bisection << ",\n\t"
                                 << "Estimate fraction: " << miser_params_.estimate_frac << ",\n\t"
                                 << "α-value: " << miser_params_.alpha << ",\n\t"
                                 << "Dither: " << miser_params_.dither << ".";

    // launch the full integration
    int res = gsl_monte_miser_integrate(function_.get(),
                                        &xlow_[0],
                                        &xhigh_[0],
                                        function_->dim,
                                        ncvg_,
                                        gsl_rng_.get(),
                                        miser_state_.get(),
                                        &result,
                                        &abserr);

    if (res != GSL_SUCCESS)
      throw CG_FATAL("Integrator:integrate") << "Error while performing the integration!\n\t"
                                             << "GSL error: " << gsl_strerror(res) << ".";

    result_ = result;
    err_result_ = abserr;
  }

  ParametersDescription IntegratorMISER::description() {
    auto desc = IntegratorGSL::description();
    desc.setDescription("MISER adaptive importance sampling integrator");
    desc.add<int>("numFunctionCalls", 50000)
        .setDescription("Number of function calls per phase space point evaluation");
    desc.add<double>("estimateFraction", 0.1);
    desc.add<int>("minCalls", 16 * 10);
    desc.add<int>("minCallsPerBisection", 32 * 16 * 10);
    desc.add<double>("alpha", 2.);
    desc.add<double>("dither", 0.1);
    return desc;
  }
}  // namespace cepgen

REGISTER_INTEGRATOR("MISER", IntegratorMISER)
