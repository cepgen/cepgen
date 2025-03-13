/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2025  Laurent Forthomme
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

#ifndef CepGen_Utils_FunctionWrapper_h
#define CepGen_Utils_FunctionWrapper_h

#include <functional>

#include "CepGen/Core/ParametersList.h"

namespace cepgen::utils {
  /// Wrapper to a 1-dimensional function with optional parameters
  class FunctionWrapper {
  public:
    explicit FunctionWrapper(const std::function<double(double)>& func) : func_(func) {}
    explicit FunctionWrapper(const std::function<double(double, const ParametersList&)>& func) : func_params_(func) {}
    explicit FunctionWrapper(const std::function<double(double, void*)>& func) : func_obj_(func) {}

    /// Call the function with a user collection of parameters
    inline double operator()(double x, const ParametersList& params = ParametersList()) const {
      if (func_params_)
        return func_params_(x, params);
      return func_(x);
    }
    /// Call the function with an unspecified object as parameters
    inline double operator()(double x, void* obj) const {
      if (func_obj_)
        return func_obj_(x, obj);
      return func_(x);
    }
    /// Call the function with a templated object as parameters
    template <typename T>
    inline double operator()(double x, const T& obj) const {
      if (func_obj_)
        return func_obj_(x, dynamic_cast<void*>(&obj));
      if (func_params_)
        return func_params_(x, ParametersList().set("parameter", obj));
      return func_(x);
    }

    inline operator const std::function<double(double)>&() { return func_; }
    inline operator const std::function<double(double, void*)>&() { return func_obj_; }

  private:
    std::function<double(double)> func_{nullptr};  ///< Reference to the parameters-less functor
    std::function<double(double, const ParametersList&)> func_params_{nullptr};  ///< Reference to the functor
    std::function<double(double, void*)> func_obj_{nullptr};                     ///< Reference to the functor
  };
}  // namespace cepgen::utils

#endif
