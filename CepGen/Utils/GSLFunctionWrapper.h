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

#ifndef CepGen_Utils_GSLFunctionWrapper_h
#define CepGen_Utils_GSLFunctionWrapper_h

#include <gsl/gsl_math.h>

#include <memory>

#include "CepGen/Utils/FunctionWrapper.h"

namespace cepgen::utils {
  /// GSL wrapper to define a functor as a GSL-digestible functional
  class GSLFunctionWrapper : public gsl_function {
  public:
    /// Utility to build a gsl_function pointer from a functional
    static std::unique_ptr<gsl_function> build(const FunctionWrapper&, void*);
    /// Utility to build a gsl_function pointer from a functional
    static std::unique_ptr<gsl_function> build(const FunctionWrapper&, const ParametersList& = ParametersList());

  private:
    explicit GSLFunctionWrapper(const FunctionWrapper&, const ParametersList&, void* j = nullptr);
    static double eval(double x, void* params);  ///< Static integrable functional
    const FunctionWrapper function_wrapper_;
    const ParametersList& params_;
    void* object_{nullptr};
  };
}  // namespace cepgen::utils

#endif
