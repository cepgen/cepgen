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

#ifndef CepGen_Utils_GSLFunctionsWrappers_h
#define CepGen_Utils_GSLFunctionsWrappers_h

#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>

#include <memory>

#include "CepGen/Utils/FunctionWrapper.h"

namespace cepgen::utils {
  /// GSL wrapper to define a functor as a GSL-digestible functional
  class GSLFunctionWrapper : public gsl_function {
  public:
    /// Utility to build a gsl_function pointer from a functional
    static std::unique_ptr<gsl_function> build(const FunctionWrapper& func, void* obj) {
      return std::unique_ptr<gsl_function>(new GSLFunctionWrapper(func, ParametersList(), obj));
    }
    /// Utility to build a gsl_function pointer from a functional
    static std::unique_ptr<gsl_function> build(const FunctionWrapper& func,
                                               const ParametersList& params = ParametersList()) {
      return std::unique_ptr<gsl_function>(new GSLFunctionWrapper(func, params, nullptr));
    }

  private:
    explicit GSLFunctionWrapper(const FunctionWrapper& func, const ParametersList& plist, void* obj = nullptr)
        : func_(func), params_(plist), obj_(obj) {
      function = &GSLFunctionWrapper::eval;
      params = this;
    }
    /// Static integrable functional
    static double eval(double x, void* params) {
      auto* wrp = static_cast<GSLFunctionWrapper*>(params);
      if (wrp->obj_)
        return wrp->func_(x, wrp->obj_);
      if (!wrp->params_.empty())
        return wrp->func_(x, wrp->params_);
      return wrp->func_(x);
    }
    const FunctionWrapper func_;
    const ParametersList& params_;
    void* obj_{nullptr};
  };

  /// GSL wrapper to define a functor as an integrable functional
  /// \tparam F functor member signature
  template <typename F>
  class GSLMonteFunctionWrapper : public gsl_monte_function {
  public:
    explicit GSLMonteFunctionWrapper(const F& func, size_t num_dimensions) : func_(func) {
      f = &GSLMonteFunctionWrapper::eval;
      dim = num_dimensions;
      params = this;
    }
    /// Utility to build a gsl_monte_function pointer from a functional and phase space size
    static std::unique_ptr<gsl_monte_function> build(const F& func, size_t num_dimensions) {
      return std::unique_ptr<gsl_monte_function>(new GSLMonteFunctionWrapper<decltype(func)>(func, num_dimensions));
    }

  private:
    /// Static integrable functional
    static double eval(double* x, size_t num_dimensions, void* params) {
      return static_cast<GSLMonteFunctionWrapper*>(params)->func_(x, num_dimensions, params);
    }
    /// Reference to the functor
    const F& func_;
  };
}  // namespace cepgen::utils

#endif
