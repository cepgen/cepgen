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

#ifndef CepGenAddOns_PythonWrapper_PythonTypes_h
#define CepGenAddOns_PythonWrapper_PythonTypes_h

#include <Python.h>

#include <memory>
#include <string>
#include <vector>

#include "CepGen/Core/ParametersList.h"
#include "CepGen/Utils/Limits.h"

#if PY_MAJOR_VERSION < 3
#define PYTHON2
#endif

namespace cepgen {
  namespace python {
    struct ObjectPtrDeleter {
      void operator()(PyObject* obj) { Py_DECREF(obj); }
    };
    typedef std::unique_ptr<PyObject, ObjectPtrDeleter> ObjectPtr;

    ObjectPtr importModule(const std::string&);

    template <typename T>
    bool is(PyObject* obj);

    template <typename T>
    T get(PyObject* obj);
    template <typename T>
    ObjectPtr set(const T&);

    template <typename T>
    bool isVector(PyObject* obj);
    template <typename T>
    std::vector<T> getVector(PyObject* obj);

    template <typename T>
    ObjectPtr newTuple(const std::vector<T>&);

    template <>
    bool is<bool>(PyObject* obj);
    template <>
    ObjectPtr set<bool>(const bool&);

    template <>
    bool is<int>(PyObject* obj);

    template <>
    bool is<long>(PyObject* obj);

    template <>
    int get<int>(PyObject* obj);
    template <>
    ObjectPtr set<int>(const int&);

    template <>
    unsigned long get<unsigned long>(PyObject* obj);

    template <>
    bool is<ParametersList>(PyObject* obj);
    template <>
    ParametersList get<ParametersList>(PyObject* obj);
    template <>
    ObjectPtr set<ParametersList>(const ParametersList&);

    template <>
    bool is<double>(PyObject* obj);
    template <>
    double get<double>(PyObject* obj);
    template <>
    ObjectPtr set<double>(const double&);

    template <>
    bool is<std::string>(PyObject* obj);
    template <>
    std::string get<std::string>(PyObject* obj);
    template <>
    ObjectPtr set<std::string>(const std::string&);

    template <>
    bool is<Limits>(PyObject* obj);
    template <>
    Limits get<Limits>(PyObject* obj);
  }  // namespace python
}  // namespace cepgen

#endif
