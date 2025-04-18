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

#include <gsl/gsl_rng.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Integration/GSLIntegrator.h"
#include "CepGen/Integration/Integrand.h"
#include "CepGen/Modules/RandomGeneratorFactory.h"
#include "CepGen/Utils/RandomGenerator.h"

using namespace cepgen;

GSLIntegrator::GSLIntegrator(const ParametersList& params)
    : Integrator(params),
      random_generator_(RandomGeneratorFactory::get().build(steer<ParametersList>("randomGenerator"))) {
  CG_DEBUG("GSLIntegrator") << "Random numbers generator: " << gsl_rng_name(random_generator_->engine<gsl_rng>())
                            << ".";
}

void GSLIntegrator::prepare(Integrand& integrand, const std::vector<Limits>& range) {
  // specify the integrand through the GSL wrapper
  function_ = [&integrand](double* x, size_t num_dimensions, void*) -> double {
    return integrand.eval(std::vector(x, x + num_dimensions));
  };
  if (gsl_function_ = utils::GSLMonteFunctionWrapper<decltype(function_)>::build(function_, integrand.size());
      !gsl_function_)
    throw CG_FATAL("GSLIntegrator:prepare") << "Integrand was not properly set.";
  if (gsl_function_->dim <= 0)
    throw CG_FATAL("GSLIntegrator:prepare") << "Invalid phase space dimension: " << gsl_function_->dim << ".";
  CG_DEBUG("GSLIntegrator:prepare") << "Number of integration dimensions: " << gsl_function_->dim << ".";

  if (const auto ndim = integrand.size(); ndim > 0) {  // set the integration range
    if (range.size() < ndim)
      throw CG_FATAL("GSLIntegrator:prepare")
          << "Insufficient number of limits (" << range << ") provided for dim-" << ndim << " integrand.";
    x_low_.resize(ndim);
    x_high_.resize(ndim);
    for (size_t i = 0; i < ndim; ++i)
      x_low_[i] = range.at(i).min(), x_high_[i] = range.at(i).max();
  }
}

ParametersDescription GSLIntegrator::description() {
  auto desc = Integrator::description();
  desc.add("randomGenerator", RandomGeneratorFactory::get().describeParameters("gsl"))
      .setDescription("type of random number generator to use for integration");
  return desc;
}
