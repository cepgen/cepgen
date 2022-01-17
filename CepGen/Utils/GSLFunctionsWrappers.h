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

namespace cepgen {
  namespace utils {
    /// GSL wrapper to define a functor as a GSL-digestible functional
    class GSLFunctionWrapper : public gsl_function {
    public:
      explicit GSLFunctionWrapper(const std::function<double(double)>& func) : func_(func) {
        function = &GSLFunctionWrapper::eval;
        params = this;
      }
      /// Utility to build a gsl_function pointer from a functional
      static std::unique_ptr<gsl_function> build(const std::function<double(double)>& func) {
        return std::unique_ptr<gsl_function>(new utils::GSLFunctionWrapper(func));
      }

    private:
      /// Static integrable functional
      static double eval(double x, void* params) { return static_cast<GSLFunctionWrapper*>(params)->func_(x); }
      /// Reference to the functor
      const std::function<double(double)>& func_;
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
