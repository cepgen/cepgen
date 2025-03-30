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

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/FunctionalFactory.h"
#include "CepGen/Utils/String.h"
#include "CepGenPython/Error.h"
#include "CepGenPython/Functional.h"

using namespace std::string_literals;

namespace cepgen::python {
  Functional::Functional(const ParametersList& params)
      : utils::Functional(params),
        environment_(new Environment(steer<ParametersList>("environment"))),
        name_(steerName()) {
    const auto cmd = "from math import *\n"s + "def " + steer<std::string>("functionName") + "("s +
                     utils::merge(vars_, ",") + ") -> float:\n" + "\treturn " +
                     utils::replaceAll(expression_, {{"^", "**"}}) + "\n";
    CG_DEBUG("python:Functional") << "Will compile Python expression:\n" << cmd;
    if (mod_ = ObjectPtr::defineModule("functional", cmd); !mod_)
      throw CG_ERROR("python:Functional") << "Failed to initialise the functional parser module.";
    try {
      func_ = mod_.attribute(steer<std::string>("functionName"));
      if (!func_ || !PyCallable_Check(func_.get()))
        throw PY_ERROR << "Failed to retrieve/cast the object to a Python functional.";
      if (const auto function_code = func_.attribute("__code__"); function_code) {
        if (const auto argument_names_attribute = function_code.attribute("co_varnames");
            argument_names_attribute && argument_names_attribute.isVector<std::string>()) {
          for (const auto& argument_name : argument_names_attribute.vector<std::string>())
            arguments_.emplace_back(argument_name);
          CG_DEBUG("python:Functional") << "List of arguments unpacked for function '" << name_ << "': " << arguments_
                                        << ".";
        } else
          CG_WARNING("python:Functional") << "Failed to retrieve argument names for function '" << name_ << "'.";
      } else
        CG_WARNING("python:Functional") << "Failed to retrieve code for function '" << name_ << "'.";
    } catch (const Error& err) {
      throw CG_ERROR("python:Functional")
          << "Failed to initialise the Python functional with \"" << expression_ << "\".\n"
          << err.message();
    }
  }

  Functional::Functional(const ObjectPtr& obj)
      : utils::Functional(ParametersList()), name_(obj.attribute("__name__").value<std::string>()), func_(obj.get()) {
    // Python environment is not needed, as it is already assumed to be present (if a Python object is given as an argument...)
    CG_DEBUG("python:Functional") << "Functional '" << name_ << "' parsed from object.";
    if (const auto code = ObjectPtr::wrap(PyFunction_GetCode(func_.get())); code) {
      CG_DEBUG("python:Functional") << "Functional has an associated code.";
      if (const auto arg_count = code.attribute("co_argcount"); arg_count && arg_count.is<int>()) {
        CG_DEBUG("python:Functional") << "Retrieved " << utils::s("argument", arg_count.value<int>(), true) << ".";
        for (int i = 0; i < arg_count.value<int>(); ++i) {
          vars_.emplace_back(utils::format("var_%d", i));
          values_.emplace_back(0.);
        }
      }
    } else
      CG_WARNING("python:Functional") << "Python code object was not retrieved from function '" << name_
                                      << "' object. Cannot count the arguments.";
  }

  double Functional::eval() const {
    const auto get_value = [this](const ObjectPtr& return_value) -> double {
      if (!return_value)
        throw PY_ERROR << "Invalid return type for function '" << name_ << "' call: " << return_value.get() << ".";
      if (return_value.is<double>())
        return return_value.value<double>();
      if (return_value.isVector<double>()) {
        if (const auto vec = return_value.vector<double>(); !vec.empty()) {
          if (vec.size() > 1)
            CG_WARNING("python:Functional") << "Invalid size for return vector of function '" << name_
                                            << "': " << vec.size() << ". Values: " << vec << ".";
          return vec.at(0);
        }
        throw PY_ERROR << "Empty result vector returned from function '" << name_ << "'.";
      }
      throw PY_ERROR << "Invalid return type for function '" << name_ << "' call: " << return_value << ".";
    };
    try {
      if (values_.size() == 1) {  // single-argument function is a bit simpler to handle
        if (const auto value = func_.call(values_.at(0)); value)
          return get_value(value);
      } else if (const auto func_arguments = ObjectPtr::tupleFromVector(values_); func_arguments) {
        if (const auto value = func_.call(func_arguments); value)
          return get_value(value);
      } else
        throw PY_ERROR << "Invalid functions argument building: " << values_ << ".";
    } catch (const Error& err) {
      throw CG_ERROR("python:Functional:eval")
          << "Failed to call the function '" << name_ << "' with arguments=" << values_ << ".\n"
          << err.message();
    }
    throw CG_ERROR("python:Functional:eval") << "Failed to build a tuple for the arguments.";
  }

  ParametersDescription Functional::description() {
    auto desc = utils::Functional::description();
    desc.setDescription("Python mathematical expression evaluator");
    desc.add("functionName", "custom_functional"s)
        .setDescription("Python function name (in case multiple instance have to be declared in a same environment)");
    return desc;
  }
}  // namespace cepgen::python
using PythonFunctional = cepgen::python::Functional;
REGISTER_FUNCTIONAL("python", PythonFunctional);
