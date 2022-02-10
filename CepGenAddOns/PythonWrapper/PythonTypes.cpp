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
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Utils/String.h"
#include "CepGenAddOns/PythonWrapper/PythonUtils.h"

#if PY_MAJOR_VERSION < 3
#define PYTHON2
#endif

namespace cepgen {
  namespace python {
    //------------------------------------------------------------------
    // typed retrieval helpers
    //------------------------------------------------------------------

    template <>
    bool is<int>(PyObject* obj) {
      if (!obj)
        error("Failed to retrieve integer object!");
#ifdef PYTHON2
      return PyInt_Check(obj);
#else
      return PyLong_Check(obj);
#endif
    }

    template <>
    bool is<bool>(PyObject* obj) {
      if (!obj)
        error("Failed to retrieve boolean object!");
      return PyBool_Check(obj);
    }

    template <>
    PyObject* set<bool>(const bool& val) {
      return PyBool_FromLong(val);
    }

    template <>
    bool is<long>(PyObject* obj) {
#ifdef PYTHON2
      return PyInt_Check(obj) || PyLong_Check(obj);
#else
      return PyLong_Check(obj);
#endif
    }

    template <>
    int get<int>(PyObject* obj) {
      if (!is<int>(obj))
        throw CG_ERROR("PythonHandler:get") << "Object has invalid type: integer \"" << obj->ob_type->tp_name << "\".";
#ifdef PYTHON2
      return PyInt_AsLong(obj);
#else
      return PyLong_AsLong(obj);
#endif
    }

    template <>
    PyObject* set<int>(const int& val) {
#ifdef PYTHON2
      return PyInt_FromLong(val);
#else
      return PyLong_FromLong(val);
#endif
    }

    template <>
    unsigned long get<unsigned long>(PyObject* obj) {
      if (!is<long>(obj))
        throw CG_ERROR("PythonHandler:get")
            << "Object has invalid type: unsigned long \"" << obj->ob_type->tp_name << "\".";
#ifdef PYTHON2
      return PyInt_AsUnsignedLongMask(obj);
#else
      if (!PyLong_Check(obj))
        error(utils::format("Object has invalid type: unsigned long != %s", obj->ob_type->tp_name));
      return PyLong_AsUnsignedLong(obj);
#endif
    }

    template <>
    long long get<long long>(PyObject* obj) {
      if (!is<long>(obj))
        throw CG_ERROR("PythonHandler:get")
            << "Object has invalid type: long long != \"" << obj->ob_type->tp_name << "\".";
      return PyLong_AsLongLong(obj);
    }

    template <>
    bool is<double>(PyObject* obj) {
      if (!obj)
        error("Failed to retrieve float object!");
      return PyFloat_Check(obj);
    }

    template <>
    double get<double>(PyObject* obj) {
      if (!is<double>(obj))
        throw CG_ERROR("PythonHandler:get")
            << "Object has invalid type: double != \"" << obj->ob_type->tp_name << "\".";
      return PyFloat_AsDouble(obj);
    }

    template <>
    PyObject* set<double>(const double& val) {
      return PyFloat_FromDouble(val);
    }

    template <>
    bool is<std::string>(PyObject* obj) {
      if (!obj)
        error("Failed to retrieve string object!");
#ifdef PYTHON2
      return PyString_Check(obj);
#else
      return PyUnicode_Check(obj);
#endif
    }

    template <>
    std::string get<std::string>(PyObject* obj) {
      if (!is<std::string>(obj))
        throw CG_ERROR("PythonHandler:get")
            << "Object has invalid type: string != \"" << obj->ob_type->tp_name << "\".";
      return decode(obj);
    }

    template <>
    PyObject* set<std::string>(const std::string& val) {
#ifdef PYTHON2
      return PyString_FromString(val.c_str());
#else
      return PyUnicode_FromString(val.c_str());
#endif
    }

    template <>
    bool is<Limits>(PyObject* obj) {
      if (!obj)
        error("Failed to retrieve limits object!");
      if (!isVector<double>(obj))
        return false;
      const size_t size = getVector<double>(obj).size();
      return size == 1 || size == 2;
    }

    template <>
    Limits get<Limits>(PyObject* obj) {
      if (!is<Limits>(obj))
        throw CG_ERROR("PythonHandler:get")
            << "Object has invalid type: limits != \"" << obj->ob_type->tp_name << "\".";
      const auto vec = getVector<double>(obj);
      if (vec.size() == 1)
        return Limits{vec.at(0)};
      return Limits{vec.at(0), vec.at(1)};
    }

    template <>
    bool is<ParametersList>(PyObject* obj) {
      if (!obj)
        error("Failed to retrieve parameters list object!");
      return PyDict_Check(obj);
    }

    template <typename T>
    bool isVector(PyObject* obj) {
      if (!obj)
        return false;
      if (!PyTuple_Check(obj) && !PyList_Check(obj))
        return false;
      const bool tuple = PyTuple_Check(obj);
      if ((tuple ? PyTuple_Size(obj) : PyList_Size(obj)) == 0)
        return true;
      auto* pfirst = tuple ? PyTuple_GetItem(obj, 0) : PyList_GetItem(obj, 0);
      if (!is<T>(pfirst))
        return false;
      return true;
    }

    template <typename T>
    std::vector<T> getVector(PyObject* obj) {
      if (!isVector<T>(obj))
        throw CG_ERROR("PythonHandler:get")
            << "Object has invalid type: list/tuple != \"" << obj->ob_type->tp_name << "\".";
      std::vector<T> vec;
      const bool tuple = PyTuple_Check(obj);
      const Py_ssize_t num_entries = tuple ? PyTuple_Size(obj) : PyList_Size(obj);
      //--- check every single element inside the list/tuple
      for (Py_ssize_t i = 0; i < num_entries; ++i) {
        auto* pit = tuple ? PyTuple_GetItem(obj, i) : PyList_GetItem(obj, i);
        if (!is<T>(pit))
          error("Mixed types detected in vector");
        vec.emplace_back(get<T>(pit));
      }
      return vec;
    }

    template <typename T>
    PyObject* newTuple(const std::vector<T>& vec) {
      auto* tuple = PyTuple_New(vec.size());
      for (size_t i = 0; i < vec.size(); ++i) {
        auto* value = set<T>(vec.at(i));
        PyTuple_SetItem(tuple, i, value);
        Py_DECREF(value);
      }
      return tuple;
    }

    template <>
    ParametersList get<ParametersList>(PyObject* obj) {
      if (!is<ParametersList>(obj))
        throw CG_ERROR("PythonHandler:get")
            << "Object has invalid type: parameters list != \"" << obj->ob_type->tp_name << "\".";
      ParametersList out;
      PyObject *pkey = nullptr, *pvalue = nullptr;
      Py_ssize_t pos = 0;
      while (PyDict_Next(obj, &pos, &pkey, &pvalue)) {
        const std::string skey = is<std::string>(pkey) ? get<std::string>(pkey)
                                 : is<int>(pkey)       ? std::to_string(get<int>(pkey))  // integer-type key
                                                       : "invalid";
        if (is<bool>(pvalue))
          out.set<bool>(skey, get<int>(pvalue));
        else if (is<int>(pvalue))
          out.set<int>(skey, get<int>(pvalue));
        else if (is<double>(pvalue))
          out.set<double>(skey, get<double>(pvalue));
        else if (is<std::string>(pvalue))
          out.set<std::string>(skey, get<std::string>(pvalue));
        else if (is<ParametersList>(pvalue))
          out.set<ParametersList>(skey, get<ParametersList>(pvalue));
        else if (PyTuple_Check(pvalue) || PyList_Check(pvalue)) {  // vector
          if (isVector<int>(pvalue))
            out.set<std::vector<int> >(skey, getVector<int>(pvalue));
          else if (isVector<double>(pvalue)) {
            if (is<Limits>(pvalue))
              out.set<Limits>(skey, get<Limits>(pvalue));
            out.set<std::vector<double> >(skey, getVector<double>(pvalue));
          } else if (isVector<std::string>(pvalue))
            out.set<std::vector<std::string> >(skey, getVector<std::string>(pvalue));
          else  //if (isVector<ParametersList>(pvalue))
            out.set<std::vector<ParametersList> >(skey, getVector<ParametersList>(pvalue));
        } else if (Py_IS_TYPE(pvalue, &_PyNone_Type)) {
          out.set<std::string>(skey, "None");
        } else {
          error(utils::format("Invalid object (%s) retrieved as parameters list value!", pvalue->ob_type->tp_name));
        }
      }
      return out;
    }

    template <>
    PyObject* set<ParametersList>(const ParametersList& plist) {
      auto* obj = PyDict_New();
      for (const auto& key : plist.keys(true)) {
        PyObject* val{nullptr};
        if (plist.has<bool>(key))
          val = set<bool>(plist.get<bool>(key));
        else if (plist.has<int>(key))
          val = set<int>(plist.get<int>(key));
        else if (plist.has<double>(key))
          val = set<double>(plist.get<double>(key));
        else if (plist.has<std::string>(key))
          val = set<std::string>(plist.get<std::string>(key));
        else if (plist.has<ParametersList>(key))
          val = set<ParametersList>(plist.get<ParametersList>(key));
        else if (plist.has<Limits>(key)) {
          const auto& lim = plist.get<Limits>(key);
          val = newTuple<double>({lim.min(), lim.max()});
        } else if (plist.has<std::vector<int> >(key))
          val = newTuple<int>(plist.get<std::vector<int> >(key));
        else if (plist.has<std::vector<double> >(key))
          val = newTuple<double>(plist.get<std::vector<double> >(key));
        else if (plist.has<std::vector<std::string> >(key))
          val = newTuple<std::string>(plist.get<std::vector<std::string> >(key));
        else
          error("Parameters list has an untranslatable object for key=" + key);
        if (val) {
          PyDict_SetItem(obj, encode(key), set<int>(plist.get<int>(key)));
          Py_DECREF(val);
        }
      }
      return obj;
    }

    template PyObject* newTuple<bool>(const std::vector<bool>&);
    template PyObject* newTuple<int>(const std::vector<int>&);
    template PyObject* newTuple<double>(const std::vector<double>&);
    template PyObject* newTuple<std::string>(const std::vector<std::string>&);
  }  // namespace python
}  // namespace cepgen
