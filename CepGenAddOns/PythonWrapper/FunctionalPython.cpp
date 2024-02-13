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

// clang-format off
#include "CepGenAddOns/PythonWrapper/Environment.h"
#include "CepGenAddOns/PythonWrapper/Error.h"
#include "CepGenAddOns/PythonWrapper/PythonUtils.h"
// clang-format on

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/FunctionalFactory.h"
#include "CepGen/Utils/Functional.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  namespace utils {
    class FunctionalPython final : public Functional {
    public:
      explicit FunctionalPython(const ParametersList&);
      double eval() const;

      static ParametersDescription description();

    private:
      python::Environment env_{ParametersList()};
      python::ObjectPtr mod_;
      python::ObjectPtr func_{nullptr};
    };

    FunctionalPython::FunctionalPython(const ParametersList& params)
        : Functional(params), mod_(PyModule_New("functional")) {
      PyModule_AddStringConstant(mod_.get(), "__file__", "");
      auto* local = PyModule_GetDict(mod_.get());  // borrowed
      python::ObjectPtr math(PyRun_String("from math import *", Py_file_input, local, local));
      auto expr = utils::replace_all(expression_, {{"^", "**"}});
      std::ostringstream os;
      os << "def custom_functional(";
      std::string sep;
      for (const auto& var : vars_)
        os << sep << var << ": float", sep = ", ";
      os << ") -> float:\n"
         << "\treturn " << expr << "\n";
      CG_DEBUG("FunctionalPython") << "Will compile Python expression:\n" << os.str();
      try {
        python::ObjectPtr value(PyRun_String(os.str().c_str(), Py_file_input, local, local));  // new
        if (!value)
          throw PY_ERROR;
        if (func_ = mod_.attribute("custom_functional"); !func_ || !PyCallable_Check(func_.get()))
          throw CG_ERROR("FunctionalPython") << "Failed to retrieve/cast the object to a Python functional.";
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
      return desc;
    }
  }  // namespace utils
}  // namespace cepgen
typedef cepgen::utils::FunctionalPython PythonFunctional;
REGISTER_FUNCTIONAL("python", PythonFunctional);
