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

#ifndef CepGen_Utils_GSLMonteFunctionWrapper_h
#define CepGen_Utils_GSLMonteFunctionWrapper_h

#include <gsl/gsl_monte.h>

#include <memory>

namespace cepgen::utils {
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
    const F& func_;  ///< Reference to the functor
  };
}  // namespace cepgen::utils

#endif
