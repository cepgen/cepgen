/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2023  Laurent Forthomme
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
#include "CepGenAddOns/PythonWrapper/Error.h"
#include "CepGenAddOns/PythonWrapper/PythonTypes.h"
// clang-format on

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  namespace python {
    void ObjectPtrDeleter::operator()(PyObject* obj) {
      CG_DEBUG("Python:ObjectPtrDeleter").log([&obj](auto& log) {
        log << "Destroying object at addr 0x" << obj << " (";
#if PY_VERSION_HEX >= 0x03110000
        auto* type = Py_TYPE(obj);
        if (type)
          log << "type: " << get<std::string>(PyType_GetName(type)) << ", ";
#endif
        log << "reference count: " << Py_REFCNT(obj) << ")";
      });
      Py_DECREF(obj);
    }

    std::ostream& operator<<(std::ostream& os, const ObjectPtr& ptr) {
      os << "PyObject{";
      auto repr = ObjectPtr(PyObject_Str(ptr.get()));  // new
      if (repr)
        os << get<std::string>(repr);
      return os << "}";
    }

    ObjectPtr importModule(const std::string& mod_name) {
      return ObjectPtr(PyImport_ImportModule(mod_name.c_str()));  // new
    }

    ObjectPtr defineModule(const std::string& mod_name, const std::string& code) {
      auto mod = ObjectPtr(PyImport_AddModule(mod_name.data()));
      if (!mod)
        throw PY_ERROR << "Failed to add the module.";
      auto* local_dict = PyModule_GetDict(mod.get());
      if (!local_dict)
        throw PY_ERROR << "Failed to retrieve the local dictionary from module.";
      auto run = ObjectPtr(PyRun_String(code.data(), Py_file_input, local_dict, local_dict));
      CG_DEBUG("Python:defineModule") << "New '" << mod_name << "' module initialised from Python code parsing.\n"
                                      << "List of attributes: " << getVector<std::string>(PyObject_Dir(mod.get()))
                                      << ".";
      return mod;
    }

    //------------------------------------------------------------------
    // typed retrieval helpers
    //------------------------------------------------------------------

    template <>
    ObjectPtr set<PyObject*>(PyObject* const& obj) {
      return ObjectPtr(obj);
    }

    template <>
    bool is<int>(PyObject* obj) {
      if (!obj)
        throw CG_ERROR("Python:is") << "Failed to retrieve integer object.";
#ifdef PYTHON2
      return PyInt_Check(obj);
#else
      return PyLong_Check(obj);
#endif
    }

    template <>
    bool is<bool>(PyObject* obj) {
      if (!obj)
        throw CG_ERROR("Python:is") << "Failed to retrieve boolean object.";
      return PyBool_Check(obj);
    }

    template <>
    bool get<bool>(PyObject* obj) {
      return PyObject_IsTrue(obj);
    }

    template <>
    ObjectPtr set<bool>(const bool& val) {
      return ObjectPtr(PyBool_FromLong(val));
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
        throw CG_ERROR("Python:get") << "Object has invalid type: integer != \"" << obj->ob_type->tp_name << "\".";
#ifdef PYTHON2
      return PyInt_AsLong(obj);
#else
      return PyLong_AsLong(obj);
#endif
    }

    template <>
    ObjectPtr set<int>(const int& val) {
#ifdef PYTHON2
      return ObjectPtr(PyInt_FromLong(val));
#else
      return ObjectPtr(PyLong_FromLong(val));
#endif
    }

    template <>
    unsigned long get<unsigned long>(PyObject* obj) {
      if (!is<long>(obj))
        throw CG_ERROR("Python:get") << "Object has invalid type: unsigned long != \"" << obj->ob_type->tp_name
                                     << "\".";
#ifdef PYTHON2
      return PyInt_AsUnsignedLongMask(obj);
#else
      return PyLong_AsUnsignedLong(obj);
#endif
    }

    template <>
    long long get<long long>(PyObject* obj) {
      if (!is<long>(obj))
        throw CG_ERROR("Python:get") << "Object has invalid type: long long != \"" << obj->ob_type->tp_name << "\".";
      return PyLong_AsLongLong(obj);
    }

    template <>
    bool is<double>(PyObject* obj) {
      if (!obj)
        throw CG_ERROR("Python:is") << "Failed to retrieve float object.";
      return PyFloat_Check(obj);
    }

    template <>
    double get<double>(PyObject* obj) {
      if (!is<double>(obj))
        throw CG_ERROR("Python:get") << "Object has invalid type: double != \"" << obj->ob_type->tp_name << "\".";
      return PyFloat_AsDouble(obj);
    }

    template <>
    ObjectPtr set<double>(const double& val) {
      return ObjectPtr(PyFloat_FromDouble(val));
    }

    template <>
    bool is<std::string>(PyObject* obj) {
      if (!obj)
        throw CG_ERROR("Python:is") << "Failed to retrieve string object.";
#ifdef PYTHON2
      return PyString_Check(obj);
#else
      return PyUnicode_Check(obj) || PyBytes_Check(obj);
#endif
    }

    template <>
    std::string get<std::string>(PyObject* obj) {
      if (!is<std::string>(obj))
        throw CG_ERROR("Python:get") << "Object has invalid type: string != \"" << obj->ob_type->tp_name << "\".";
#ifdef PYTHON2
      return PyString_AsString(obj);
#else
      if (PyUnicode_Check(obj))
        return PyUnicode_AsUTF8(obj);
      else  // if (PyBytes_Check(obj))
        return strdup(PyBytes_AS_STRING(obj));
#endif
    }

    template <>
    ObjectPtr set<std::string>(const std::string& val) {
#ifdef PYTHON2
      return ObjectPtr(PyString_FromString(val.c_str()));
#else
      return ObjectPtr(PyUnicode_FromString(val.c_str()));
#endif
    }

    template <>
    bool is<Limits>(PyObject* obj) {
      if (!obj)
        throw CG_ERROR("Python:is") << "Failed to retrieve limits object.";
      if (!isVector<double>(obj))
        return false;
      const size_t size = getVector<double>(obj).size();
      return size == 1 || size == 2;
    }

    template <>
    Limits get<Limits>(PyObject* obj) {
      if (!is<Limits>(obj))
        throw CG_ERROR("Python:get") << "Object has invalid type: limits != \"" << obj->ob_type->tp_name << "\".";
      const auto vec = getVector<double>(obj);
      if (vec.size() == 1)
        return Limits{vec.at(0)};
      return Limits{vec.at(0), vec.at(1)};
    }

    template <>
    ObjectPtr set<Limits>(const Limits& val) {
      return newTuple(std::vector<double>{val.min(), val.max()});
    }

    template <>
    bool is<ParametersList>(PyObject* obj) {
      if (!obj)
        throw CG_ERROR("Python:is") << "Failed to retrieve parameters list object.";
      return PyDict_Check(obj);
    }

    template <typename T>
    bool isVector(PyObject* obj) {
      if (!obj)
        throw CG_ERROR("Python:isVector") << "Object is not a vector ; in fact, it is not even an object.";
      const bool tuple = PyTuple_Check(obj), list = PyList_Check(obj);
      if (!tuple && !list)
        return false;
      const auto size = tuple ? PyTuple_Size(obj) : list ? PyList_Size(obj) : 0;
      if (size == 0)
        return true;
      auto* pfirst = tuple ? PyTuple_GetItem(obj, 0) : list ? PyList_GetItem(obj, 0) : nullptr;
      if (!pfirst)
        return false;
      if (!is<T>(pfirst)) {  // only allow same-type tuples/lists
        CG_WARNING("Python::isVector") << "Wrong object type unpacked from tuple/list: (python)"
                                       << pfirst->ob_type->tp_name << " != (c++)" << utils::demangle(typeid(T).name())
                                       << ".";
        return false;
      }
      return true;
    }

    template <typename T>
    std::vector<T> getVector(PyObject* obj) {
      if (!isVector<T>(obj))
        throw CG_ERROR("Python:getVector")
            << "Object has invalid type: list/tuple != \"" << obj->ob_type->tp_name << "\".";
      std::vector<T> vec;
      const bool tuple = PyTuple_Check(obj);
      const Py_ssize_t num_entries = tuple ? PyTuple_Size(obj) : PyList_Size(obj);
      //--- check every single element inside the list/tuple
      for (Py_ssize_t i = 0; i < num_entries; ++i) {
        auto* pit = tuple ? PyTuple_GetItem(obj, i) : PyList_GetItem(obj, i);
        if (!is<T>(pit))
          throw CG_ERROR("Python:getVector") << "Mixed types detected in vector.";
        vec.emplace_back(get<T>(pit));
      }
      return vec;
    }

    template <typename T>
    ObjectPtr newTuple(const std::vector<T>& vec) {
      ObjectPtr tuple(PyTuple_New(vec.size()));
      for (size_t i = 0; i < vec.size(); ++i)
        PyTuple_SetItem(tuple.get(), i, set<T>(vec.at(i)).release());
      return tuple;
    }

    template <>
    ObjectPtr newTuple(const std::vector<PyObject*>& vec) {
      ObjectPtr tuple(PyTuple_New(vec.size()));
      for (size_t i = 0; i < vec.size(); ++i)
        PyTuple_SetItem(tuple.get(), i, vec.at(i));
      return tuple;
    }

    template <>
    ParametersList get<ParametersList>(PyObject* obj) {
      if (!is<ParametersList>(obj))
        throw CG_ERROR("Python:get") << "Object has invalid type: parameters list != \"" << obj->ob_type->tp_name
                                     << "\".";
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
        } else if (pvalue == Py_None) {
          out.set<std::string>(skey, "None");
        } else {
          CG_WARNING("PythonTypes") << "Invalid object (" << pvalue->ob_type->tp_name << ") retrieved for key=" << skey
                                    << " when unpacking a dictionary/parameters list.";
        }
      }
      return out;
    }

    template <>
    ObjectPtr set<ParametersList>(const ParametersList& plist) {
      ObjectPtr obj(PyDict_New());
      for (const auto& key : plist.keys(true)) {
        if (plist.has<bool>(key))
          PyDict_SetItem(obj.get(), set(key).release(), set(plist.get<bool>(key)).release());
        else if (plist.has<int>(key))
          PyDict_SetItem(obj.get(), set(key).release(), set(plist.get<int>(key)).release());
        else if (plist.has<double>(key))
          PyDict_SetItem(obj.get(), set(key).release(), set(plist.get<double>(key)).release());
        else if (plist.has<std::string>(key))
          PyDict_SetItem(obj.get(), set(key).release(), set(plist.get<std::string>(key)).release());
        else if (plist.has<ParametersList>(key))
          PyDict_SetItem(obj.get(), set(key).release(), set(plist.get<ParametersList>(key)).release());
        else if (plist.has<Limits>(key)) {
          const auto& lim = plist.get<Limits>(key);
          PyDict_SetItem(obj.get(), set(key).release(), newTuple<double>({lim.min(), lim.max()}).release());
        } else if (plist.has<std::vector<int> >(key))
          PyDict_SetItem(obj.get(), set(key).release(), newTuple<int>(plist.get<std::vector<int> >(key)).release());
        else if (plist.has<std::vector<double> >(key))
          PyDict_SetItem(
              obj.get(), set(key).release(), newTuple<double>(plist.get<std::vector<double> >(key)).release());
        else if (plist.has<std::vector<std::string> >(key))
          PyDict_SetItem(obj.get(),
                         set(key).release(),
                         newTuple<std::string>(plist.get<std::vector<std::string> >(key)).release());
        else
          throw CG_ERROR("Python:set") << "Parameters list has an untranslatable object for key=" << key;
      }
      return obj;
    }

    template bool isVector<ParametersList>(PyObject*);
    template ObjectPtr newTuple<bool>(const std::vector<bool>&);
    template ObjectPtr newTuple<int>(const std::vector<int>&);
    template ObjectPtr newTuple<double>(const std::vector<double>&);
    template ObjectPtr newTuple<std::string>(const std::vector<std::string>&);
  }  // namespace python
}  // namespace cepgen
