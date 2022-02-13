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

#ifndef CepGenAddOns_PythonWrapper_PythonUtils_h
#define CepGenAddOns_PythonWrapper_PythonUtils_h

#include <Python.h>

#include <string>
#include <vector>

#include "CepGen/Core/ParametersList.h"
#include "CepGen/Utils/Limits.h"

#define PY_ERROR cepgen::python::Error(__FUNC__)

namespace cepgen {
  namespace python {
    class Error final : public LoggedException {
    public:
      explicit Error(const std::string&);
      ~Error();

    private:
      PyObject* ptype_{nullptr};
      PyObject* pvalue_{nullptr};
      PyObject* ptraceback_obj_{nullptr};
    };
    struct PyObject_deleter {
      void operator()(PyObject* obj) { Py_DECREF(obj); }
    };
    typedef std::unique_ptr<PyObject, PyObject_deleter> ObjectPtr;

    std::string pythonPath(const std::string&);
    PyObject* element(PyObject*, const std::string&);
    ObjectPtr encode(const std::string&);
    std::string decode(PyObject* obj);

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

    void fillParameter(PyObject* parent, const char* key, bool& out);
    void fillParameter(PyObject* parent, const char* key, int& out);
    void fillParameter(PyObject* parent, const char* key, unsigned long& out);
    void fillParameter(PyObject* parent, const char* key, unsigned int& out);
    void fillParameter(PyObject* parent, const char* key, double& out);
    void fillParameter(PyObject* parent, const char* key, std::string& out);
    void fillParameter(PyObject* parent, const char* key, Limits& out);
    void fillParameter(PyObject* parent, const char* key, std::vector<int>& out);
    void fillParameter(PyObject* parent, const char* key, std::vector<double>& out);
    void fillParameter(PyObject* parent, const char* key, std::vector<std::string>& out);
    void fillParameter(PyObject* parent, const char* key, ParametersList& out);
    void fillParameter(PyObject* parent, const char* key, std::vector<ParametersList>& out);

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
