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

#include "CepGen/Utils/GSLFunctionWrapper.h"

using namespace cepgen::utils;

GSLFunctionWrapper::GSLFunctionWrapper(const FunctionWrapper& function_wrapper,
                                       const ParametersList& parameters,
                                       void* object)
    : function_wrapper_(function_wrapper), params_(parameters), object_(object) {
  function = &GSLFunctionWrapper::eval;
  params = this;
}

std::unique_ptr<gsl_function> GSLFunctionWrapper::build(const FunctionWrapper& function_wrapper, void* object) {
  return std::unique_ptr<gsl_function>(new GSLFunctionWrapper(function_wrapper, ParametersList(), object));
}

std::unique_ptr<gsl_function> GSLFunctionWrapper::build(const FunctionWrapper& function_wrapper,
                                                        const ParametersList& parameters) {
  return std::unique_ptr<gsl_function>(new GSLFunctionWrapper(function_wrapper, parameters, nullptr));
}

double GSLFunctionWrapper::eval(double x, void* params) {
  auto* wrapper = reinterpret_cast<GSLFunctionWrapper*>(params);
  if (wrapper->object_)
    return wrapper->function_wrapper_(x, wrapper->object_);
  if (!wrapper->params_.empty())
    return wrapper->function_wrapper_(x, wrapper->params_);
  return wrapper->function_wrapper_(x);
}
