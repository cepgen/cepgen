/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2018-2024  Laurent Forthomme
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

#include "CepGen/Core/ParametersList.h"
#include "CepGen/Utils/Limits.h"
#include "CepGen/Utils/Message.h"
#include "CepGen/Utils/String.h"
#include "CepGenPython/Error.h"
#include "CepGenPython/ObjectPtr.h"

using namespace std::string_literals;

namespace cepgen::python {
  void obj_deleter(PyObject* obj) {
    CG_DEBUG("python:ObjectPtrDeleter").log([&obj](auto& log) {
      log << "Destroying object at addr 0x" << obj << " (";
#if PY_VERSION_HEX >= 0x03110000
      if (auto* type = Py_TYPE(obj); type)
        log << "type: " << value<std::string>(PyType_GetName(type)) << ", ";
#endif
      log << "reference count: " << Py_REFCNT(obj) << ")";
    });
    Py_DECREF(obj);
  }

  ObjectPtr::ObjectPtr(PyObject* obj, bool wrap_only)
      : PyObjectPtr(obj, wrap_only ? [](PyObject*) { /* do not dereference if only wrapping */ } : obj_deleter) {}

  ObjectPtr ObjectPtr::wrap(PyObject* obj) {
    ObjectPtr ptr(obj, true);
    return ptr;
  }

  // generic, un-specialised makers/getters

  template <typename T>
  ObjectPtr ObjectPtr::make(const T&) {
    throw CG_FATAL("ObjectPtr:make") << "Type specialisation is not implemented.";
  }

  template <typename T>
  bool ObjectPtr::is() const {
    throw CG_FATAL("ObjectPtr:is") << "Type specialisation is not implemented.";
  }

  template <typename T>
  T ObjectPtr::value() const {
    throw CG_FATAL("ObjectPtr:value") << "Type specialisation is not implemented.";
  }

  // type specialisations

  template <>
  bool ObjectPtr::is<int>() const {
    CG_ASSERT(get());
#ifdef PYTHON2
    return PyInt_Check(get());
#else
    return PyLong_Check(get());
#endif
  }

  template <>
  bool ObjectPtr::is<bool>() const {
    CG_ASSERT(get());
    return PyBool_Check(get());
  }

  template <>
  bool ObjectPtr::is<long>() const {
    CG_ASSERT(get());
#ifdef PYTHON2
    return PyInt_Check(get()) || PyLong_Check(get());
#else
    return PyLong_Check(get());
#endif
  }

  template <>
  bool ObjectPtr::is<double>() const {
    CG_ASSERT(get());
    return PyFloat_Check(get());
  }

  template <>
  bool ObjectPtr::is<std::string>() const {
    CG_ASSERT(get());
#ifdef PYTHON2
    return PyString_Check(get());
#else
    return PyUnicode_Check(get()) || PyBytes_Check(get());
#endif
  }

  template <>
  bool ObjectPtr::is<Limits>() const {
    if (!isVector<double>())
      return false;
    if (const auto size = vector<double>().size(); size == 1 || size == 2)
      return true;
    return false;
  }

  template <>
  bool ObjectPtr::is<ParametersList>() const {
    CG_ASSERT(get());
    return PyDict_Check(get());
  }

  template <typename T>
  bool ObjectPtr::isVector() const {
    CG_ASSERT(get());
    const bool tuple = PyTuple_Check(get()), list = PyList_Check(get());
    if (!tuple && !list)  // only accept 'tuples' and 'lists'
      return false;
    if (const auto size = tuple ? PyTuple_Size(get()) : list ? PyList_Size(get()) : 0; size == 0)
      return true;
    const auto first = ObjectPtr::wrap(tuple  ? PyTuple_GetItem(get(), 0) /* borrowed */
                                       : list ? PyList_GetItem(get(), 0)  /* borrowed */
                                              : nullptr);
    if (!first)
      return false;
    if (!first.is<T>()) {  // only allow same-type tuples/lists
      CG_DEBUG("python:ObjectPtr:isVector")
          << "Wrong object type unpacked from tuple/list: (python)" << first->ob_type->tp_name << " != (c++)"
          << utils::demangle(typeid(T).name()) << ".";
      return false;
    }
    return true;
  }

  template <typename T>
  std::vector<T> ObjectPtr::vector() const {
    if (!isVector<T>())
      throw CG_ERROR("python::ObjectPtr:vector")
          << "Object has invalid type: list/tuple != \"" << get()->ob_type->tp_name << "\".";
    std::vector<T> vec;
    const bool tuple = PyTuple_Check(get());
    const Py_ssize_t num_entries = tuple ? PyTuple_Size(get()) : PyList_Size(get());
    //--- check every single element inside the list/tuple
    for (Py_ssize_t i = 0; i < num_entries; ++i) {
      auto pit =
          ObjectPtr::wrap(tuple ? PyTuple_GetItem(get(), i) /* borrowed */ : PyList_GetItem(get(), i) /* borrowed */);
      if (!pit.is<T>())
        throw CG_ERROR("python::ObjectPtr:vector") << "Mixed types detected in vector.";
      vec.emplace_back(pit.value<T>());
    }
    return vec;
  }

  template <>
  bool ObjectPtr::value<bool>() const {
    CG_ASSERT(get());
    return PyObject_IsTrue(get());
  }

  template <>
  int ObjectPtr::value<int>() const {
    if (!is<int>())
      throw CG_ERROR("Python:get") << "Object has invalid type: integer != \"" << get()->ob_type->tp_name << "\".";
#ifdef PYTHON2
    return PyInt_AsLong(get());
#else
    return PyLong_AsLong(get());
#endif
  }

  template <>
  unsigned long ObjectPtr::value<unsigned long>() const {
    if (!is<long>())
      throw CG_ERROR("Python:get") << "Object has invalid type: unsigned long != \"" << get()->ob_type->tp_name
                                   << "\".";
#ifdef PYTHON2
    return PyInt_AsUnsignedLongMask(get());
#else
    return PyLong_AsUnsignedLong(get());
#endif
  }

  template <>
  long long ObjectPtr::value<long long>() const {
    if (!is<long>())
      throw CG_ERROR("Python:get") << "Object has invalid type: long long != \"" << get()->ob_type->tp_name << "\".";
    return PyLong_AsLongLong(get());
  }

  template <>
  double ObjectPtr::value<double>() const {
    if (!is<double>())
      throw CG_ERROR("Python:get") << "Object has invalid type: double != \"" << get()->ob_type->tp_name << "\".";
    return PyFloat_AsDouble(get());
  }

  template <>
  std::string ObjectPtr::value<std::string>() const {
    if (!is<std::string>())
      throw CG_ERROR("Python:get") << "Object has invalid type: string != \"" << get()->ob_type->tp_name << "\".";
#ifdef PYTHON2
    return PyString_AsString(get());
#else
    if (PyUnicode_Check(get()))
      return PyUnicode_AsUTF8(get());
    if (auto* str_buf = ::strdup(PyBytes_AS_STRING(get())); str_buf) {
      auto out = std::string(str_buf);
      free(str_buf);
      return out;
    }
    throw CG_ERROR("Python:get") << "Failed to retrieve a string buffer from object.";
#endif
  }

  template <>
  Limits ObjectPtr::value<Limits>() const {
    if (!is<Limits>())
      throw CG_ERROR("Python:get") << "Object has invalid type: limits != \"" << get()->ob_type->tp_name << "\".";
    const auto vec = vector<double>();
    if (vec.size() == 1)
      return Limits{vec.at(0)};
    return Limits{vec.at(0), vec.at(1)};
  }

  template <>
  ParametersList ObjectPtr::value<ParametersList>() const {
    if (!is<ParametersList>())
      throw CG_ERROR("Python:get") << "Object has invalid type: parameters list != \"" << get()->ob_type->tp_name
                                   << "\".";
    ParametersList out;
    Py_ssize_t pos = 0;
    PyObject *pkey{nullptr}, *pvalue{nullptr};
    while (PyDict_Next(get(), &pos, &pkey, &pvalue)) {
      const auto key = ObjectPtr::wrap(pkey), val = ObjectPtr::wrap(pvalue);
      const std::string skey = key.is<std::string>() ? key.value<std::string>()
                               : key.is<int>()       ? std::to_string(key.value<int>())  // integer-type key
                                                     : "invalid";
      if (val.is<bool>())
        out.set(skey, static_cast<bool>(val.value<int>()));
      else if (val.is<int>())
        out.set(skey, val.value<int>());
      else if (val.is<double>())
        out.set(skey, val.value<double>());
      else if (val.is<std::string>())
        out.set(skey, val.value<std::string>());
      else if (val.is<ParametersList>())
        out.set(skey, val.value<ParametersList>());
      else if (PyTuple_Check(pvalue) || PyList_Check(pvalue)) {  // vector
        if (val.isVector<int>())
          out.set(skey, val.vector<int>());
        else if (val.isVector<double>()) {
          if (val.is<Limits>())
            out.set(skey, val.value<Limits>());
          out.set(skey, val.vector<double>());
        } else if (val.isVector<std::string>())
          out.set(skey, val.vector<std::string>());
        else if (val.isVector<Limits>())
          out.set(skey, val.vector<Limits>());
        else  //if (val.isVector<ParametersList>())
          out.set(skey, val.vector<ParametersList>());
      } else if (pvalue == Py_None) {
        out.set(skey, "None"s);
      } else {
        CG_WARNING("PythonTypes") << "Invalid object (" << pvalue->ob_type->tp_name << ") retrieved for key=" << skey
                                  << " when unpacking a dictionary/parameters list.";
      }
    }
    return out;
  }

  template <>
  ObjectPtr ObjectPtr::make<PyObject*>(PyObject* const& obj) {
    return ObjectPtr(obj);
  }

  template <>
  ObjectPtr ObjectPtr::make<int>(const int& val) {
#ifdef PYTHON2
    return ObjectPtr(PyInt_FromLong(val));
#else
    return ObjectPtr(PyLong_FromLong(val));
#endif
  }

  template <>
  ObjectPtr ObjectPtr::make<bool>(const bool& val) {
    return ObjectPtr(PyBool_FromLong(val));
  }

  template <>
  ObjectPtr ObjectPtr::make<double>(const double& val) {
    return ObjectPtr(PyFloat_FromDouble(val));
  }

  template <>
  ObjectPtr ObjectPtr::make<std::string>(const std::string& val) {
#ifdef PYTHON2
    return ObjectPtr(PyString_FromString(val.c_str()));
#else
    return ObjectPtr(PyUnicode_FromString(val.c_str()));
#endif
  }

  template <>
  ObjectPtr ObjectPtr::make<Limits>(const Limits& val) {
    return ObjectPtr::tupleFromVector(std::vector<double>{val.min(), val.max()});
  }

  template <>
  ObjectPtr ObjectPtr::make<ParametersList>(const ParametersList& plist) {
    ObjectPtr obj(PyDict_New());
    for (const auto& key : plist.keys(true)) {
      if (plist.has<bool>(key))
        PyDict_SetItem(obj.get(), make(key).release(), make(plist.get<bool>(key)).release());
      else if (plist.has<int>(key))
        PyDict_SetItem(obj.get(), make(key).release(), make(plist.get<int>(key)).release());
      else if (plist.has<double>(key))
        PyDict_SetItem(obj.get(), make(key).release(), make(plist.get<double>(key)).release());
      else if (plist.has<std::string>(key))
        PyDict_SetItem(obj.get(), make(key).release(), make(plist.get<std::string>(key)).release());
      else if (plist.has<ParametersList>(key))
        PyDict_SetItem(obj.get(), make(key).release(), make(plist.get<ParametersList>(key)).release());
      else if (plist.has<Limits>(key)) {
        const auto& lim = plist.get<Limits>(key);
        PyDict_SetItem(
            obj.get(), make(key).release(), ObjectPtr::tupleFromVector<double>({lim.min(), lim.max()}).release());
      } else if (plist.has<std::vector<int> >(key))
        PyDict_SetItem(obj.get(),
                       make(key).release(),
                       ObjectPtr::tupleFromVector<int>(plist.get<std::vector<int> >(key)).release());
      else if (plist.has<std::vector<double> >(key))
        PyDict_SetItem(obj.get(),
                       make(key).release(),
                       ObjectPtr::tupleFromVector<double>(plist.get<std::vector<double> >(key)).release());
      else if (plist.has<std::vector<std::string> >(key))
        PyDict_SetItem(obj.get(),
                       make(key).release(),
                       ObjectPtr::tupleFromVector<std::string>(plist.get<std::vector<std::string> >(key)).release());
      else
        throw PY_ERROR << "Parameters list has an untranslatable object for key=" << key;
    }
    return obj;
  }

  template <typename T>
  ObjectPtr ObjectPtr::tupleFromVector(const std::vector<T>& vec) {
    ObjectPtr tuple(PyTuple_New(vec.size()));
    for (size_t i = 0; i < vec.size(); ++i)
      PyTuple_SetItem(tuple.get(), i, make<T>(vec.at(i)).release());
    return tuple;
  }

  template <>
  ObjectPtr ObjectPtr::tupleFromVector(const std::vector<PyObject*>& vec) {
    ObjectPtr tuple(PyTuple_New(vec.size()));
    for (size_t i = 0; i < vec.size(); ++i)
      PyTuple_SetItem(tuple.get(), i, vec.at(i));
    return tuple;
  }

  template ObjectPtr ObjectPtr::tupleFromVector<bool>(const std::vector<bool>&);
  template ObjectPtr ObjectPtr::tupleFromVector<int>(const std::vector<int>&);
  template ObjectPtr ObjectPtr::tupleFromVector<double>(const std::vector<double>&);
  template ObjectPtr ObjectPtr::tupleFromVector<std::string>(const std::vector<std::string>&);
  template ObjectPtr ObjectPtr::tupleFromVector<Limits>(const std::vector<Limits>&);

  std::ostream& operator<<(std::ostream& os, const ObjectPtr& ptr) {
    os << "PyObject{";
    if (auto repr = ObjectPtr(PyObject_Str(ptr.get())); repr)  // new
      os << repr.value<std::string>();
    return os << "}";
  }

  ObjectPtr ObjectPtr::call(const ObjectPtr& args) const {
    return ObjectPtr(PyObject_CallObject(get(), args.get()) /* new */);
  }

  ObjectPtr ObjectPtr::attribute(const std::string& attr) const {
    if (PyObject_HasAttrString(get(), attr.c_str()) != 1)
      return ObjectPtr(nullptr);
    return ObjectPtr(PyObject_GetAttrString(get(), attr.c_str()));  // new
  }

  ObjectPtr ObjectPtr::importModule(const std::string& mod_name) {
    return ObjectPtr(PyImport_Import(ObjectPtr::make<std::string>(mod_name).get()));  // new
  }

  ObjectPtr ObjectPtr::defineModule(const std::string& mod_name, const std::string& code) {
    auto mod = ObjectPtr(PyImport_AddModule(mod_name.data()));
    if (!mod)
      throw PY_ERROR << "Failed to add the module.";
    auto* local_dict = PyModule_GetDict(mod.get());
    if (!local_dict)
      throw PY_ERROR << "Failed to retrieve the local dictionary from module.";
    ObjectPtr::wrap(PyRun_String(code.data(), Py_file_input, local_dict, local_dict));
    CG_DEBUG("Python:defineModule") << "New '" << mod_name << "' module initialised from Python code parsing.\n"
                                    << "List of attributes: "
                                    << ObjectPtr(PyObject_Dir(mod.get())).vector<std::string>() << ".";
    return mod;
  }
}  // namespace cepgen::python
