/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2022  Laurent Forthomme
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

#ifndef CepGen_Utils_GSLFunctionsWrappers_h
#define CepGen_Utils_GSLFunctionsWrappers_h

#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>

#include <memory>

#include "CepGen/Core/ParametersList.h"

namespace cepgen {
  namespace utils {
    /// Wrapper to a 1-dimensional function with optional parameters
    class Function1D {
    public:
      Function1D(const std::function<double(double)>& func) : func_(func), func_params_(nullptr) {}
      Function1D(const std::function<double(double, const ParametersList&)>& func)
          : func_(nullptr), func_params_(func) {}
      double operator()(double x, const ParametersList& params = ParametersList()) const {
        if (func_params_)
          return func_params_(x, params);
        return func_(x);
      }

    private:
      /// Reference to the parameters-less functor
      std::function<double(double)> func_;
      /// Reference to the functor
      std::function<double(double, const ParametersList&)> func_params_;
    };

    /// GSL wrapper to define a functor as a GSL-digestible functional
    class GSLFunctionWrapper : public gsl_function {
    public:
      /// Utility to build a gsl_function pointer from a functional
      static std::unique_ptr<gsl_function> build(const Function1D& func,
                                                 const ParametersList& params = ParametersList()) {
        return std::unique_ptr<gsl_function>(new utils::GSLFunctionWrapper(func, params));
      }

    private:
      explicit GSLFunctionWrapper(const Function1D& func, const ParametersList& plist) : func_(func), params_(plist) {
        function = &GSLFunctionWrapper::eval;
        params = this;
      }
      /// Static integrable functional
      static double eval(double x, void* params) {
        auto* wrp = static_cast<GSLFunctionWrapper*>(params);
        return wrp->func_(x, wrp->params_);
      }
      const Function1D func_;
      const ParametersList& params_;
    };

    /// GSL wrapper to define a functor as an integrable functional
    /// \tparam F functor member signature
    template <typename F>
    class GSLMonteFunctionWrapper : public gsl_monte_function {
    public:
      explicit GSLMonteFunctionWrapper(const F& func, size_t ndim) : func_(func) {
        f = &GSLMonteFunctionWrapper::eval;
        dim = ndim;
        params = this;
      }
      /// Utility to build a gsl_monte_function pointer from a functional and phase space size
      static std::unique_ptr<gsl_monte_function> build(const F& func, size_t ndim) {
        return std::unique_ptr<gsl_monte_function>(new utils::GSLMonteFunctionWrapper<decltype(func)>(func, ndim));
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
