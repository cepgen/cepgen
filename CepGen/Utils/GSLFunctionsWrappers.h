#ifndef CepGen_Utils_GSLFunctionsWrappers_h
#define CepGen_Utils_GSLFunctionsWrappers_h

#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>

namespace cepgen {
  namespace utils {
    /**
     * A GSL wrapper to define a functor as a GSL-digestible functional
     * \tparam F functor member signature
     */
    template <typename F>
    class GSLFunctionWrapper : public gsl_function {
    public:
      GSLFunctionWrapper(const F& func) : func_(func) {
        function = &GSLFunctionWrapper::eval;
        params = this;
      }

    private:
      /// Static integrable functional
      static double eval(double* x, void* params) { return static_cast<GSLFunctionWrapper*>(params)->func_(x, params); }
      /// Reference to the functor
      const F& func_;
    };
    /**
     * A GSL wrapper to define a functor as an integrable functional
     * \tparam F functor member signature
     */
    template <typename F>
    class GSLMonteFunctionWrapper : public gsl_monte_function {
    public:
      GSLMonteFunctionWrapper(const F& func, size_t ndim) : func_(func) {
        f = &GSLMonteFunctionWrapper::eval;
        dim = ndim;
        params = this;
      }

    private:
      /// Static integrable functional
      static double eval(double* x, size_t ndim, void* params) {
        return static_cast<GSLMonteFunctionWrapper*>(params)->func_(x, ndim, params);
      }
      /// Reference to the functor
      const F& func_;
    };
  }  // namespace utils
}  // namespace cepgen

#endif
