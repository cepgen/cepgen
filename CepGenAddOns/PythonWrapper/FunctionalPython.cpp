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

#include <Python.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/FunctionalFactory.h"
#include "CepGen/Utils/Functional.h"
#include "CepGenAddOns/PythonWrapper/PythonUtils.h"

namespace cepgen {
  namespace utils {
    class FunctionalPython final : public Functional {
    public:
      explicit FunctionalPython(const ParametersList&);
      ~FunctionalPython();

      static ParametersDescription description();

      double eval(const std::vector<double>&) const;

    private:
      mutable PyObject* value_;
      PyObject* func_;
      PyObject* mod_;
    };

    FunctionalPython::FunctionalPython(const ParametersList& params) : Functional(params) {
      Py_Initialize();
      auto* global = PyDict_New();
      mod_ = PyModule_New("functional");
      PyModule_AddStringConstant(mod_, "__file__", "");
      auto* local = PyModule_GetDict(mod_);
      std::ostringstream os;
      os << "from math import *\n"
         << "def custom_functional(";
      std::string sep;
      for (const auto& var : vars_)
        os << sep << var, sep = ", ";
      os << "):\n"
         << "\treturn " << expression_ << "\n";
      CG_LOG << os.str();
      value_ = PyRun_String(os.str().c_str(), Py_file_input, global, local);
      if (!value_)
        python::error("Failed to initialise the Python functional with \"" + expression_ + "\".");
      Py_DECREF(value_);
      func_ = PyObject_GetAttrString(mod_, "custom_functional");
      if (!func_ || !PyCallable_Check(func_))
        python::error("Failed to retrieve/cast the object to a Python functional.");
    }

    FunctionalPython::~FunctionalPython() {
      Py_DECREF(value_);
      Py_XDECREF(func_);
      Py_DECREF(mod_);
      Py_Finalize();
    }

    double FunctionalPython::eval(const std::vector<double>& x) const {
      auto* args = PyTuple_New(x.size());
      for (size_t i = 0; i < x.size(); ++i) {
        value_ = PyFloat_FromDouble(x.at(i));
        PyTuple_SetItem(args, i, value_);
      }
      value_ = PyObject_CallObject(func_, args);
      Py_DECREF(args);
      //if (!PyFloat_Check(value_))
      //  throw CG_ERROR("FunctionalPython") << "Return object is not a float.";
      auto ret = PyFloat_AsDouble(value_);
      CG_LOG << "ret:" << ret;
      return ret;
    }

    ParametersDescription FunctionalPython::description() {
      auto desc = Functional::description();
      desc.setDescription("Python mathematical expression evaluator");
      return desc;
    }
  }  // namespace utils
}  // namespace cepgen

REGISTER_FUNCTIONAL("python", FunctionalPython)
