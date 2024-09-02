/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2024  Laurent Forthomme
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

namespace cepgen {
  GSLIntegrator::GSLIntegrator(const ParametersList& params) : Integrator(params) {
    CG_DEBUG("Integrator:build") << "Random numbers generator: " << gsl_rng_name(rnd_gen_->engine<gsl_rng>()) << ".";
  }

  void GSLIntegrator::setIntegrand(Integrand& integrand) {
    //--- specify the integrand through the GSL wrapper
    function_ = [&integrand](double* x, size_t num_dimensions, void*) -> double {
      return integrand.eval(std::vector<double>(x, x + num_dimensions));
    };
    gsl_function_ = utils::GSLMonteFunctionWrapper<decltype(function_)>::build(function_, integrand.size());
    if (!gsl_function_)
      throw CG_FATAL("GSLIntegrator:setIntegrand") << "Integrand was not properly set.";
    if (gsl_function_->dim <= 0)
      throw CG_FATAL("GSLIntegrator:setIntegrand") << "Invalid phase space dimension: " << gsl_function_->dim << ".";

    CG_DEBUG("GSLIntegrator:setIntegrand") << "Number of integration dimensions: " << gsl_function_->dim << ".";

    checkLimits(integrand);  // check the integration bounds
  }

  void GSLIntegrator::setLimits(const std::vector<Limits>& limits) {
    Integrator::setLimits(limits);
    x_low_.clear();
    x_high_.clear();
    for (const auto& lim : limits_) {
      x_low_.emplace_back(lim.min());
      x_high_.emplace_back(lim.max());
    }
  }

  ParametersDescription GSLIntegrator::description() {
    auto desc = Integrator::description();
    desc.add<ParametersDescription>("randomGenerator", RandomGeneratorFactory::get().describeParameters("gsl"));
    return desc;
  }
}  // namespace cepgen
