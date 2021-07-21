#include "CepGen/Integration/IntegratorGSL.h"
#include "CepGen/Integration/Integrand.h"

#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/Exception.h"

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
}  // namespace cepgen
