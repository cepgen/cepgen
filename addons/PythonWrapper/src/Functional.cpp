/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2024  Laurent Forthomme
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
      : utils::Functional(params), environment_(new Environment(steer<ParametersList>("environment"))) {
    const auto cmd = "from math import *\n"s + "def " + steer<std::string>("functionName") + "("s +
                     utils::merge(vars_, ",") + ") -> float:\n" + "\treturn " +
                     utils::replaceAll(expression_, {{"^", "**"}}) + "\n";
    CG_DEBUG("python:Functional") << "Will compile Python expression:\n" << cmd;
    mod_ = ObjectPtr::defineModule("functional", cmd);
    try {
      if (func_ = mod_.attribute(steer<std::string>("functionName")); !func_ || !PyCallable_Check(func_.get()))
        throw PY_ERROR << "Failed to retrieve/cast the object to a Python functional.";
    } catch (const Error& err) {
      throw CG_ERROR("python:Functional")
          << "Failed to initialise the Python functional with \"" << expression_ << "\".\n"
          << err.message();
    }
  }

  Functional::Functional(const ObjectPtr& obj) : utils::Functional(ParametersList()), func_(obj.get()) {
    // Python environment is not needed, as it is already assumed to be present (if a Python object is given as an argument...)
    CG_DEBUG("python:Functional") << "Functional '" << obj.attribute("__name__").value<std::string>()
                                  << "' parsed from object.";
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
      CG_WARNING("python:Functional")
          << "Python code object was not retrieved from function object. Cannot count the arguments.";
  }

  double Functional::eval() const {
    if (values_.size() == 1)
      if (const auto value = func_.call(values_.at(0)); value)
        return value.value<double>();
    if (const auto func_arguments = ObjectPtr::tupleFromVector(values_); func_arguments)
      try {
        if (const auto value = func_.call(func_arguments); value)
          return value.value<double>();
        throw PY_ERROR;
      } catch (const Error& err) {
        throw CG_ERROR("python:Functional:eval")
            << "Failed to call the function with arguments=" << func_arguments.vector<double>() << ".\n"
            << err.message();
      }
    throw CG_ERROR("python:Functional:eval") << "Failed to build a tuple for the arguments.";
  }

  ParametersDescription Functional::description() {
    auto desc = utils::Functional::description();
    desc.setDescription("Python mathematical expression evaluator");
    desc.add<std::string>("functionName", "custom_functional")
        .setDescription("Python function name (in case multiple instance have to be declared in a same environment)");
    return desc;
  }
}  // namespace cepgen::python
using PythonFunctional = cepgen::python::Functional;
REGISTER_FUNCTIONAL("python", PythonFunctional);
