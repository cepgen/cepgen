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

#include <gsl/gsl_monte_vegas.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Integration/GSLIntegrator.h"
#include "CepGen/Integration/Integrand.h"
#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGen/Utils/ProcessVariablesAnalyser.h"
#include "CepGen/Utils/String.h"

using namespace cepgen;

/// Vegas integration algorithm developed by P. Lepage, as documented in \cite Lepage:1977sw
class VegasIntegrator final : public GSLIntegrator {
public:
  explicit VegasIntegrator(const ParametersList& params)
      : GSLIntegrator(params),
        num_function_calls_(steer<int>("numFunctionCalls")),
        chi_square_cut_(steer<double>("chiSqCut")),
        treat_(steer<bool>("treat")) {
    verbosity_ = steer<int>("verbose");  // supersede the parent default verbosity level
  }

  static ParametersDescription description() {
    auto desc = GSLIntegrator::description();
    desc.setDescription("Vegas stratified sampling integrator");
    desc.add<int>("numFunctionCalls", 100'000);
    desc.add<double>("chiSqCut", 1.5);
    desc.add<bool>("treat", true).setDescription("Phase space treatment");
    desc.add<int>("iterations", 10);
    desc.add<double>("alpha", 1.25);
    desc.addAs<int, Mode>("mode", Mode::stratified);
    desc.add<std::string>("loggingOutput", "cerr");
    desc.add<int>("verbose", -1);
    return desc;
  }

  Value integrate(Integrand&) override;

  enum class Mode { importance = 1, importanceOnly = 0, stratified = -1 };
  friend std::ostream& operator<<(std::ostream& os, const Mode& mode) {
    switch (mode) {
      case Mode::importance:
        return os << "importance";
      case Mode::importanceOnly:
        return os << "importance-only";
      case Mode::stratified:
        return os << "stratified";
    }
    return os;
  }

private:
  void warmup(size_t);

  double COORD(size_t i, size_t j) const { return vegas_state_->xi[i * vegas_state_->dim + j]; }

  double eval(Integrand& integrand, const std::vector<double>& x) const override {
    if (!treat_)  // by default, no grid treatment
      return integrand.eval(x);
    //--- treatment of the integration grid
    if (r_boxes_ == 0) {
      r_boxes_ = static_cast<size_t>(std::pow(vegas_state_->bins, integrand.size()));
      x_new_.resize(integrand.size());
    }
    double w = r_boxes_;
    for (size_t j = 0; j < integrand.size(); ++j) {
      //--- find surrounding coordinates and interpolate
      const double z = x.at(j) * vegas_state_->bins;
      const auto id = static_cast<size_t>(z);  // coordinate of point before
      const double rel_pos = z - id;           // position between coordinates (norm.)
      const double bin_width = (id == 0) ? COORD(1, j) : COORD(id + 1, j) - COORD(id, j);
      //--- build new coordinate from linear interpolation
      x_new_[j] = COORD(id + 1, j) - bin_width * (1. - rel_pos);
      w *= bin_width;
    }
    return w * integrand.eval(x_new_);
  }

  const int num_function_calls_;
  const double chi_square_cut_;
  const bool treat_;  ///< Is the integrand to be smoothed for events generation?
  gsl_monte_vegas_params vegas_params_;

  /// A trivial deleter for the Vegas integrator
  struct gsl_monte_vegas_deleter {
    inline void operator()(gsl_monte_vegas_state* state) const { gsl_monte_vegas_free(state); }
  };

  /// A Vegas integrator state for integration (optional) and/or
  /// "treated" event generation
  std::unique_ptr<gsl_monte_vegas_state, gsl_monte_vegas_deleter> vegas_state_;
  mutable unsigned long long r_boxes_{0ull};
  mutable std::vector<double> x_new_;
};

Value VegasIntegrator::integrate(Integrand& integrand) {
  setIntegrand(integrand);

  //--- start by preparing the grid/state
  vegas_state_.reset(gsl_monte_vegas_alloc(gsl_function_->dim));
  gsl_monte_vegas_params_get(vegas_state_.get(), &vegas_params_);
  vegas_params_.iterations = steer<int>("iterations");
  vegas_params_.alpha = steer<double>("alpha");
  vegas_params_.verbose = verbosity_;
  vegas_params_.mode = steer<int>("mode");
  //--- output logging
  const auto& log = steer<std::string>("loggingOutput");
  if (log == "cerr")  // redirect all debugging information to the error stream
    vegas_params_.ostream = stderr;
  else if (log == "cout")  // redirect all debugging information to the standard stream
    vegas_params_.ostream = stdout;
  else
    vegas_params_.ostream = fopen(log.c_str(), "w");
  gsl_monte_vegas_params_set(vegas_state_.get(), &vegas_params_);

  //--- a bit of printout for debugging
  CG_DEBUG("Integrator:build") << "Vegas parameters:\n\t"
                               << "Number of iterations in Vegas: " << vegas_params_.iterations << ",\n\t"
                               << "α-value: " << vegas_params_.alpha << ",\n\t"
                               << "Verbosity: " << vegas_params_.verbose << ",\n\t"
                               << "Grid interpolation mode: " << static_cast<VegasIntegrator::Mode>(vegas_params_.mode)
                               << ".";
  if (!vegas_state_)
    throw CG_FATAL("Integrator:integrate") << "Vegas state not initialised!";

  //--- launch integration

  warmup(25'000);  // warmup (prepare the grid)

  // integration phase
  unsigned short chi_square = 0;
  double result, absolute_error;
  do {
    if (int res = gsl_monte_vegas_integrate(gsl_function_.get(),
                                            &x_low_[0],
                                            &x_high_[0],
                                            gsl_function_->dim,
                                            0.2 * num_function_calls_,
                                            rnd_gen_->engine<gsl_rng>(),
                                            vegas_state_.get(),
                                            &result,
                                            &absolute_error);
        res != GSL_SUCCESS)
      throw CG_FATAL("Integrator:integrate")
          << "Error at iteration #" << chi_square << " while performing the integration!\n\t"
          << "GSL error: " << gsl_strerror(res) << ".";
    CG_LOG << "\t>> at call " << (++chi_square) << ": "
           << utils::format(
                  "average = %10.6f   "
                  "sigma = %10.6f   chi2 = %4.3f.",
                  result,
                  absolute_error,
                  gsl_monte_vegas_chisq(vegas_state_.get()));
  } while (std::fabs(gsl_monte_vegas_chisq(vegas_state_.get()) - 1.) > chi_square_cut_ - 1.);
  CG_DEBUG("Integrator:integrate") << "Vegas grid information:\n\t"
                                   << "ran for " << vegas_state_->dim << " dimensions, "
                                   << "and generated " << vegas_state_->bins_max << " bins.\n\t"
                                   << "Integration volume: " << vegas_state_->vol << ".";

  return Value{result, absolute_error};
}

void VegasIntegrator::warmup(size_t num_calls) {
  if (!vegas_state_)
    throw CG_FATAL("Integrator:warmup") << "Vegas state not initialised!";
  // perform a first integration to warm up the grid
  double result = 0., absolute_error = 0.;
  // ensure the operation was successful
  if (int res = gsl_monte_vegas_integrate(gsl_function_.get(),
                                          &x_low_[0],
                                          &x_high_[0],
                                          gsl_function_->dim,
                                          num_calls,
                                          rnd_gen_->engine<gsl_rng>(),
                                          vegas_state_.get(),
                                          &result,
                                          &absolute_error);

      res != GSL_SUCCESS)
    throw CG_ERROR("VegasIntegrator:warmup") << "Failed to warm-up the Vegas grid.\n\t"
                                             << "GSL error: " << gsl_strerror(res) << ".";

  CG_INFO("VegasIntegrator:warmup") << "Finished the Vegas warm-up.";
}
REGISTER_INTEGRATOR("Vegas", VegasIntegrator);
