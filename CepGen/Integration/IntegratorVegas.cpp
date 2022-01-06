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

#include <gsl/gsl_monte_vegas.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Integration/GridParameters.h"
#include "CepGen/Integration/Integrand.h"
#include "CepGen/Integration/IntegratorGSL.h"
#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGen/Parameters.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  /// Vegas integration algorithm developed by P. Lepage, as documented in \cite Lepage:1977sw
  class IntegratorVegas final : public IntegratorGSL {
  public:
    explicit IntegratorVegas(const ParametersList&);

    static ParametersDescription description();

    void integrate(double&, double&) override;

    enum class Mode { importance = 1, importanceOnly = 0, stratified = -1 };

  private:
    void warmup(std::vector<double>&, std::vector<double>&, unsigned int);

    double eval(const std::vector<double>&) const override;

    const int ncvg_;
    const double chisq_cut_;
    const bool treat_;  ///< Is the integrand to be smoothed for events generation?
    gsl_monte_vegas_params vegas_params_;

    /// A trivial deleter for the Vegas integrator
    struct gsl_monte_vegas_deleter {
      inline void operator()(gsl_monte_vegas_state* state) { gsl_monte_vegas_free(state); }
    };

    /// A Vegas integrator state for integration (optional) and/or
    /// "treated" event generation
    std::unique_ptr<gsl_monte_vegas_state, gsl_monte_vegas_deleter> vegas_state_;
    mutable unsigned long long r_boxes_{0ull};
    mutable std::vector<double> x_new_;
  };

  std::ostream& operator<<(std::ostream&, const IntegratorVegas::Mode&);

  IntegratorVegas::IntegratorVegas(const ParametersList& params)
      : IntegratorGSL(params),
        ncvg_(steer<int>("numFunctionCalls")),
        chisq_cut_(steer<double>("chiSqCut")),
        treat_(steer<bool>("treat")) {
    verbosity_ = steer<int>("verbose");  // supersede the parent default verbosity level
  }

  void IntegratorVegas::integrate(double& result, double& abserr) {
    if (!initialised_) {
      //--- start by preparing the grid/state
      vegas_state_.reset(gsl_monte_vegas_alloc(function_->dim));
      gsl_monte_vegas_params_get(vegas_state_.get(), &vegas_params_);
      vegas_params_.iterations = steer<int>("iterations");
      vegas_params_.alpha = steer<double>("alpha");
      vegas_params_.verbose = verbosity_;
      vegas_params_.mode = steer<int>("mode");
      //--- output logging
      const auto& log = steer<std::string>("loggingOutput");
      if (log == "cerr")
        // redirect all debugging information to the error stream
        vegas_params_.ostream = stderr;
      else if (log == "cout")
        // redirect all debugging information to the standard stream
        vegas_params_.ostream = stdout;
      else
        vegas_params_.ostream = fopen(log.c_str(), "w");
      gsl_monte_vegas_params_set(vegas_state_.get(), &vegas_params_);

      //--- a bit of printout for debugging
      CG_DEBUG("Integrator:build") << "Vegas parameters:\n\t"
                                   << "Number of iterations in Vegas: " << vegas_params_.iterations << ",\n\t"
                                   << "α-value: " << vegas_params_.alpha << ",\n\t"
                                   << "Verbosity: " << vegas_params_.verbose << ",\n\t"
                                   << "Grid interpolation mode: " << (IntegratorVegas::Mode)vegas_params_.mode << ".";
      initialised_ = true;
    }
    if (!vegas_state_)
      throw CG_FATAL("Integrator:integrate") << "Vegas state not initialised!";

    //--- integration bounds
    std::vector<double> x_low(function_->dim, 0.), x_up(function_->dim, 1.);

    //--- launch integration

    //----- warmup (prepare the grid)
    warmup(x_low, x_up, 25000);

    //----- integration
    unsigned short it_chisq = 0;
    do {
      int res = gsl_monte_vegas_integrate(function_.get(),
                                          &x_low[0],
                                          &x_up[0],
                                          function_->dim,
                                          0.2 * ncvg_,
                                          gsl_rng_.get(),
                                          vegas_state_.get(),
                                          &result,
                                          &abserr);
      CG_LOG << "\t>> at call " << (++it_chisq) << ": "
             << utils::format(
                    "average = %10.6f   "
                    "sigma = %10.6f   chi2 = %4.3f.",
                    result,
                    abserr,
                    gsl_monte_vegas_chisq(vegas_state_.get()));
      if (res != GSL_SUCCESS)
        throw CG_FATAL("Integrator:integrate")
            << "Error at iteration #" << it_chisq << " while performing the integration!\n\t"
            << "GSL error: " << gsl_strerror(res) << ".";
    } while (fabs(gsl_monte_vegas_chisq(vegas_state_.get()) - 1.) > chisq_cut_ - 1.);
    CG_DEBUG("Integrator:integrate") << "Vegas grid information:\n\t"
                                     << "ran for " << vegas_state_->dim << " dimensions, "
                                     << "and generated " << vegas_state_->bins_max << " bins.\n\t"
                                     << "Integration volume: " << vegas_state_->vol << ".";

    result_ = result;
    err_result_ = abserr;
  }

  void IntegratorVegas::warmup(std::vector<double>& x_low, std::vector<double>& x_up, unsigned int ncall) {
    if (!vegas_state_)
      throw CG_FATAL("Integrator:warmup") << "Vegas state not initialised!";

    // perform a first integration to warm up the grid
    double result = 0., abserr = 0.;
    int res = gsl_monte_vegas_integrate(function_.get(),
                                        &x_low[0],
                                        &x_up[0],
                                        function_->dim,
                                        ncall,
                                        gsl_rng_.get(),
                                        vegas_state_.get(),
                                        &result,
                                        &abserr);

    // ensure the operation was successful
    if (res != GSL_SUCCESS)
      throw CG_ERROR("Integrator:vegas") << "Failed to warm-up the Vegas grid.\n\t"
                                         << "GSL error: " << gsl_strerror(res) << ".";

    CG_INFO("Integrator:vegas") << "Finished the Vegas warm-up.";
  }

  double IntegratorVegas::eval(const std::vector<double>& x) const {
    if (!integrand_)
      throw CG_FATAL("Integrator:vegas") << "Invalid integrand specified!";
    //--- by default, no grid treatment
    if (!treat_)
      return integrand_->eval(x);
    //--- treatment of the integration grid
    if (r_boxes_ == 0) {
      r_boxes_ = (size_t)std::pow(vegas_state_->bins, integrand_->size());
      x_new_.resize(integrand_->size());
    }
    auto COORD = [](gsl_monte_vegas_state* s, size_t i, size_t j) { return s->xi[i * s->dim + j]; };
    double w = r_boxes_;
    for (size_t j = 0; j < integrand_->size(); ++j) {
      //--- find surrounding coordinates and interpolate
      const double z = x.at(j) * vegas_state_->bins;
      const auto id = (size_t)z;      // coordinate of point before
      const double rel_pos = z - id;  // position between coordinates (norm.)
      const double bin_width = (id == 0) ? COORD(vegas_state_.get(), 1, j)
                                         : COORD(vegas_state_.get(), id + 1, j) - COORD(vegas_state_.get(), id, j);
      //--- build new coordinate from linear interpolation
      x_new_[j] = COORD(vegas_state_.get(), id + 1, j) - bin_width * (1. - rel_pos);
      w *= bin_width;
    }
    return w * integrand_->eval(x_new_);
  }

  std::ostream& operator<<(std::ostream& os, const IntegratorVegas::Mode& mode) {
    switch (mode) {
      case IntegratorVegas::Mode::importance:
        return os << "importance";
      case IntegratorVegas::Mode::importanceOnly:
        return os << "importance-only";
      case IntegratorVegas::Mode::stratified:
        return os << "stratified";
    }
    return os;
  }

  ParametersDescription IntegratorVegas::description() {
    auto desc = IntegratorGSL::description();
    desc.setDescription("Vegas stratified sampling integrator");
    desc.add<int>("numFunctionCalls", 50000);
    desc.add<double>("chiSqCut", 1.5);
    desc.add<bool>("treat", true).setDescription("Phase space treatment");
    desc.add<int>("iterations", 10);
    desc.add<double>("alpha", 1.5);
    desc.add<int>("mode", (int)Mode::importance);
    desc.add<std::string>("loggingOutput", "cerr");
    desc.add<int>("verbose", -1);
    return desc;
  }
}  // namespace cepgen

REGISTER_INTEGRATOR("Vegas", IntegratorVegas)
