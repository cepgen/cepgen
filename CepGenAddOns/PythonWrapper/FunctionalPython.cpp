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
      PyObject* func_;
      PyObject* mod_;
    };

    FunctionalPython::FunctionalPython(const ParametersList& params) : Functional(params) {
      Py_Initialize();
      auto* global = PyEval_GetGlobals();  // borrowed
      //auto* math = PyImport_ImportModule("math");
      //PyDict_Merge(global, PyModule_GetDict(math), true);
      mod_ = PyModule_New("functional");  // new
      PyModule_AddStringConstant(mod_, "__file__", "");
      auto* local = PyModule_GetDict(mod_);  // borrowed
      std::ostringstream os;
      os  //<< "from math import *\n"
          << "def custom_functional(";
      std::string sep;
      for (const auto& var : vars_)
        os << sep << var << ": float", sep = ", ";
      os << ") -> float:\n"
         << "\treturn " << expression_ << "\n";
      CG_DEBUG("FunctionalPython") << "Will compile Python expression:\n" << os.str();
      {
        auto* value = PyRun_String(os.str().c_str(), Py_file_input, global, local);  // new
        //PyRun_SimpleString(os.str().c_str());
        if (!value)
          python::error("Failed to initialise the Python functional with \"" + expression_ + "\".");
        CG_DEBUG("FunctionalPython") << "Python expression compilation output: " << python::get<ParametersList>(value);
        Py_DECREF(value);
      }
      func_ = PyObject_GetAttrString(mod_, "custom_functional");
      if (!func_ || !PyCallable_Check(func_))
        python::error("Failed to retrieve/cast the object to a Python functional.");
    }

    FunctionalPython::~FunctionalPython() {
      Py_XDECREF(func_);
      Py_DECREF(mod_);
      Py_Finalize();
    }

    double FunctionalPython::eval(const std::vector<double>& x) const {
      auto* args = python::newTuple(x);
      auto* value = PyObject_CallObject(func_, args);
      Py_DECREF(args);
      auto ret = python::get<double>(value);
      CG_LOG << "ret:" << ret;
      Py_DECREF(value);
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
