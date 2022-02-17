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

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/FunctionalFactory.h"
#include "CepGen/Utils/Functional.h"
#include "CepGenAddOns/PythonWrapper/PythonError.h"
#include "CepGenAddOns/PythonWrapper/PythonTypes.h"
#include "CepGenAddOns/PythonWrapper/PythonUtils.h"

namespace cepgen {
  namespace utils {
    class FunctionalPython final : public Functional {
    public:
      explicit FunctionalPython(const ParametersList&);
      double eval() const;

      static ParametersDescription description();

    private:
      python::Environment env_;
      python::ObjectPtr mod_;
      python::ObjectPtr func_;
    };

    FunctionalPython::FunctionalPython(const ParametersList& params)
        : Functional(params), mod_(PyModule_New("functional")) {
      PyModule_AddStringConstant(mod_.get(), "__file__", "");
      auto* local = PyModule_GetDict(mod_.get());  // borrowed
      python::ObjectPtr math(PyRun_String("from math import *", Py_file_input, local, local));
      std::ostringstream os;
      os << "def custom_functional(";
      std::string sep;
      for (const auto& var : vars_)
        os << sep << var << ": float", sep = ", ";
      os << ") -> float:\n"
         << "\treturn " << expression_ << "\n";
      CG_DEBUG("FunctionalPython") << "Will compile Python expression:\n" << os.str();
      {
        python::ObjectPtr value(PyRun_String(os.str().c_str(), Py_file_input, local, local));  // new
        if (!value)
          throw PY_ERROR << "Failed to initialise the Python functional with \"" << expression_ << "\".";
      }
      func_ = python::getAttribute(mod_.get(), "custom_functional");
      if (!func_ || !PyCallable_Check(func_.get()))
        throw PY_ERROR << "Failed to retrieve/cast the object to a Python functional.";
    }

    double FunctionalPython::eval() const {
      auto args = python::newTuple(values_);
      python::ObjectPtr value(PyObject_CallObject(func_.get(), args.get()));  // new
      if (!value)
        throw CG_ERROR("FunctionalPython:eval")
            << "Failed to call the function with arguments=" << python::getVector<double>(args) << ".";
      return python::get<double>(value.get());
    }

    ParametersDescription FunctionalPython::description() {
      auto desc = Functional::description();
      desc.setDescription("Python mathematical expression evaluator");
      return desc;
    }
  }  // namespace utils
}  // namespace cepgen

REGISTER_FUNCTIONAL("python", FunctionalPython)
