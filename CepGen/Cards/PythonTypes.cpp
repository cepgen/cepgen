#include "CepGen/Cards/PythonHandler.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Utils/String.h"

#if PY_MAJOR_VERSION < 3
#define PYTHON2
#endif

namespace cepgen {
  namespace card {
    //------------------------------------------------------------------
    // typed retrieval helpers
    //------------------------------------------------------------------

    template <>
    bool PythonHandler::is<int>(PyObject* obj) const {
      if (!obj)
        throwPythonError("Failed to retrieve integer object!");
#ifdef PYTHON2
      return PyInt_Check(obj);
#else
      return PyLong_Check(obj);
#endif
    }

    template <>
    bool PythonHandler::is<bool>(PyObject* obj) const {
      if (!obj)
        throwPythonError("Failed to retrieve boolean object!");
      return PyBool_Check(obj);
    }

    template <>
    bool PythonHandler::is<long>(PyObject* obj) const {
#ifdef PYTHON2
      return PyInt_Check(obj) || PyLong_Check(obj);
#else
      return PyLong_Check(obj);
#endif
    }

    template <>
    int PythonHandler::get<int>(PyObject* obj) const {
      if (!is<int>(obj))
        throw CG_ERROR("PythonHandler:get") << "Object has invalid type: integer \"" << obj->ob_type->tp_name << "\".";
#ifdef PYTHON2
      return PyInt_AsLong(obj);
#else
      return PyLong_AsLong(obj);
#endif
    }

    template <>
    unsigned long PythonHandler::get<unsigned long>(PyObject* obj) const {
      if (!is<long>(obj))
        throw CG_ERROR("PythonHandler:get")
            << "Object has invalid type: unsigned long \"" << obj->ob_type->tp_name << "\".";
#ifdef PYTHON2
      return PyInt_AsUnsignedLongMask(obj);
#else
      if (!PyLong_Check(obj))
        throwPythonError(utils::format("Object has invalid type: unsigned long != %s", obj->ob_type->tp_name));
      return PyLong_AsUnsignedLong(obj);
#endif
    }

    template <>
    long long PythonHandler::get<long long>(PyObject* obj) const {
      if (!is<long>(obj))
        throw CG_ERROR("PythonHandler:get")
            << "Object has invalid type: long long != \"" << obj->ob_type->tp_name << "\".";
      return PyLong_AsLongLong(obj);
    }

    template <>
    bool PythonHandler::is<double>(PyObject* obj) const {
      if (!obj)
        throwPythonError("Failed to retrieve float object!");
      return PyFloat_Check(obj);
    }

    template <>
    double PythonHandler::get<double>(PyObject* obj) const {
      if (!is<double>(obj))
        throw CG_ERROR("PythonHandler:get")
            << "Object has invalid type: double != \"" << obj->ob_type->tp_name << "\".";
      return PyFloat_AsDouble(obj);
    }

    template <>
    bool PythonHandler::is<std::string>(PyObject* obj) const {
      if (!obj)
        throwPythonError("Failed to retrieve string object!");
#ifdef PYTHON2
      return PyString_Check(obj);
#else
      return PyUnicode_Check(obj);
#endif
    }

    template <>
    std::string PythonHandler::get<std::string>(PyObject* obj) const {
      if (!is<std::string>(obj))
        throw CG_ERROR("PythonHandler:get")
            << "Object has invalid type: string != \"" << obj->ob_type->tp_name << "\".";
      return decode(obj);
    }

    template <>
    bool PythonHandler::is<Limits>(PyObject* obj) const {
      if (!obj)
        throwPythonError("Failed to retrieve limits object!");
      if (!isVector<double>(obj))
        return false;
      const size_t size = getVector<double>(obj).size();
      return (size == 1 || size == 2);
    }

    template <>
    Limits PythonHandler::get<Limits>(PyObject* obj) const {
      if (!is<Limits>(obj))
        throw CG_ERROR("PythonHandler:get")
            << "Object has invalid type: limits != \"" << obj->ob_type->tp_name << "\".";
      const auto vec = getVector<double>(obj);
      if (vec.size() == 1)
        return Limits{vec.at(0)};
      return Limits{vec.at(0), vec.at(1)};
    }

    template <>
    bool PythonHandler::is<ParametersList>(PyObject* obj) const {
      if (!obj)
        throwPythonError("Failed to retrieve parameters list object!");
      return PyDict_Check(obj);
    }

    template <typename T>
    bool PythonHandler::isVector(PyObject* obj) const {
      if (!obj)
        return false;
      if (!PyTuple_Check(obj) && !PyList_Check(obj))
        return false;
      const bool tuple = PyTuple_Check(obj);
      if ((tuple ? PyTuple_Size(obj) : PyList_Size(obj)) == 0)
        return true;
      PyObject* pfirst = tuple ? PyTuple_GetItem(obj, 0) : PyList_GetItem(obj, 0);
      if (!is<T>(pfirst))
        return false;
      return true;
    }

    template <typename T>
    std::vector<T> PythonHandler::getVector(PyObject* obj) const {
      if (!isVector<T>(obj))
        throw CG_ERROR("PythonHandler:get")
            << "Object has invalid type: list/tuple != \"" << obj->ob_type->tp_name << "\".";
      std::vector<T> vec;
      const bool tuple = PyTuple_Check(obj);
      const Py_ssize_t num_entries = tuple ? PyTuple_Size(obj) : PyList_Size(obj);
      //--- check every single element inside the list/tuple
      for (Py_ssize_t i = 0; i < num_entries; ++i) {
        PyObject* pit = tuple ? PyTuple_GetItem(obj, i) : PyList_GetItem(obj, i);
        if (!is<T>(pit))
          throwPythonError("Mixed types detected in vector");
        vec.emplace_back(get<T>(pit));
      }
      return vec;
    }

    template <>
    ParametersList PythonHandler::get<ParametersList>(PyObject* obj) const {
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
        if (is<int>(pvalue))
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
            out.set<std::vector<double> >(skey, getVector<double>(pvalue));
            if (is<Limits>(pvalue))
              out.set<Limits>(skey, get<Limits>(pvalue));
          } else if (isVector<std::string>(pvalue))
            out.set<std::vector<std::string> >(skey, getVector<std::string>(pvalue));
          else  //if (isVector<ParametersList>(pvalue))
            out.set<std::vector<ParametersList> >(skey, getVector<ParametersList>(pvalue));
        } else
          throwPythonError("Invalid object retrieved as parameters list value!");
      }
      return out;
    }
  }  // namespace card
}  // namespace cepgen
