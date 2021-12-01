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

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Integration/Integrand.h"
#include "CepGen/Integration/IntegratorGSL.h"
#include "CepGen/Parameters.h"
#include "CepGen/Processes/Process.h"

namespace cepgen {
  IntegratorGSL::IntegratorGSL(const ParametersList& params)
      : Integrator(params), funct_([=](double* x, size_t ndim, void*) -> double {
          return integrand_->eval(std::vector<double>(x, x + ndim));
        }) {
    //--- initialise the random number generator
    const auto& rng_type = params.get<int>("rngEngine");
    gsl_rng_type* rng_engine = nullptr;
    switch (rng_type) {
      case 0:
      default:
        rng_engine = (gsl_rng_type*)gsl_rng_mt19937;
        break;
      case 1:
        rng_engine = (gsl_rng_type*)gsl_rng_taus2;
        break;
      case 2:
        rng_engine = (gsl_rng_type*)gsl_rng_gfsr4;
        break;
      case 3:
        rng_engine = (gsl_rng_type*)gsl_rng_ranlxs0;
        break;
    }
    if (!rng_engine)
      throw CG_FATAL("Integrator:build") << "Random number generator engine not set!";

    gsl_rng_.reset(gsl_rng_alloc(rng_engine));
    gsl_rng_set(gsl_rng_.get(), seed_);

    //--- a bit of printout for debugging

    CG_DEBUG("Integrator:build") << "Random numbers generator: " << gsl_rng_name(gsl_rng_.get()) << ".\n\t"
                                 << "Seed: " << seed_ << ".";
  }

  void IntegratorGSL::setIntegrand(Integrand& integr) {
    integrand_ = &integr;
    //--- specify the integrand through the GSL wrapper
    function_.reset(new gsl_monte_function_wrapper<decltype(funct_)>(funct_, integrand_->size()));

    CG_DEBUG("Integrator:integrand") << "Number of integration dimensions: " << function_->dim << ".";

    //--- force the reinitialisation
    initialised_ = false;
  }

  double IntegratorGSL::uniform() const {
    if (!gsl_rng_)
      throw CG_FATAL("Integrator:uniform") << "Random number generator has not been initialised!";
    return gsl_rng_uniform(gsl_rng_.get());
  }

  ParametersDescription IntegratorGSL::parametersDescription() {
    auto desc = Integrator::parametersDescription();
    desc.add<int>("rngEngine", 0)
        .setDescription("Random number generator engine (0 = MT19937, 1 = Taus2, 2 = Gfsr4, 3 = RanLXS0)");
    return desc;
  }
}  // namespace cepgen
