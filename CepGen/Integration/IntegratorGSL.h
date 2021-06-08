#ifndef CepGen_Integration_IntegratorGSL_h
#define CepGen_Integration_IntegratorGSL_h

#include "CepGen/Integration/Integrator.h"

#include <memory>
#include <functional>

#include <gsl/gsl_monte.h>
#include <gsl/gsl_rng.h>

namespace cepgen {
  class IntegratorGSL : public Integrator {
  public:
    explicit IntegratorGSL(const ParametersList&);
    double uniform() const override;
    void setIntegrand(Integrand& integr) override;

  protected:
    /// A functor wrapping GSL's function footprint
    std::function<double(double*, size_t, void*)> funct_;
    /// GSL structure storing the function to be integrated by this
    /// integrator instance (along with its parameters)
    std::unique_ptr<gsl_monte_function> function_;
    /// A deleter object for GSL's random number generator
    struct gsl_rng_deleter {
      /// Destructor method for the random number generator service
      inline void operator()(gsl_rng* rng) { gsl_rng_free(rng); }
    };
    /// Instance of random number generator service
    std::unique_ptr<gsl_rng, gsl_rng_deleter> gsl_rng_;

  private:
    /**
       * \brief A GSL wrapper to define a functor as an integrable functional
       * \tparam F Class member signature
       */
    template <typename F>
    class gsl_monte_function_wrapper : public gsl_monte_function {
    public:
      gsl_monte_function_wrapper(const F& func, size_t ndim) : func_(func) {
        f = &gsl_monte_function_wrapper::eval;
        dim = ndim;
        params = this;
      }

    private:
      /// Static integrable functional
      static double eval(double* x, size_t ndim, void* params) {
        return static_cast<gsl_monte_function_wrapper*>(params)->func_(x, ndim, params);
      }
      /// Reference to the functor
      const F& func_;
    };
  };
}  // namespace cepgen

#endif
