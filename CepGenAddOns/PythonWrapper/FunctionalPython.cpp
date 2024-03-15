/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2023  Laurent Forthomme
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
#include "CepGen/Utils/Functional.h"
#include "CepGen/Utils/String.h"
#include "CepGenAddOns/PythonWrapper/Environment.h"
#include "CepGenAddOns/PythonWrapper/Error.h"
#include "CepGenAddOns/PythonWrapper/PythonUtils.h"

using namespace std::string_literals;

namespace cepgen {
  namespace utils {
    class FunctionalPython final : public Functional {
    public:
      explicit FunctionalPython(const ParametersList&);
      double eval() const;

      static ParametersDescription description();

    private:
      python::Environment env_{ParametersList()};
      python::ObjectPtr mod_{nullptr};
      python::ObjectPtr func_{nullptr};
    };

    FunctionalPython::FunctionalPython(const ParametersList& params) : Functional(params) {
      const auto cmd = "from math import *\ndef " + steer<std::string>("functionName") + "("s +
                       utils::merge(vars_, ",") + ") -> float:\n\treturn " +
                       utils::replaceAll(expression_, {{"^", "**"}}) + "\n";
      CG_DEBUG("FunctionalPython") << "Will compile Python expression:\n" << cmd;
      mod_ = python::ObjectPtr::defineModule("functional", cmd);
      try {
        if (func_ = mod_.attribute(steer<std::string>("functionName")); !func_ || !PyCallable_Check(func_.get()))
          throw PY_ERROR << "Failed to retrieve/cast the object to a Python functional.";
      } catch (const python::Error& err) {
        throw CG_ERROR("FunctionalPython")
            << "Failed to initialise the Python functional with \"" << expression_ << "\".\n"
            << err.message();
      }
    }

    double FunctionalPython::eval() const {
      auto args = python::ObjectPtr::tupleFromVector(values_);
      try {
        if (auto value = func_.call(args); value)
          return value.value<double>();
        throw PY_ERROR;
      } catch (const python::Error& err) {
        throw CG_ERROR("FunctionalPython:eval")
            << "Failed to call the function with arguments=" << args.vector<double>() << ".\n"
            << err.message();
      }
    }

    ParametersDescription FunctionalPython::description() {
      auto desc = Functional::description();
      desc.setDescription("Python mathematical expression evaluator");
      desc.add<std::string>("functionName", "custom_functional")
          .setDescription("Python function name (in case multiple instance have to be declared in a same environment)");
      return desc;
    }
  }  // namespace utils
}  // namespace cepgen
typedef cepgen::utils::FunctionalPython PythonFunctional;
REGISTER_FUNCTIONAL("python", PythonFunctional);
