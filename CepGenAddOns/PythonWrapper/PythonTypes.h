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

#ifndef CepGenAddOns_PythonWrapper_PythonTypes_h
#define CepGenAddOns_PythonWrapper_PythonTypes_h

#include <tuple>
#include <vector>

#include "CepGen/Core/ParametersList.h"
#include "CepGen/Utils/Limits.h"
#include "CepGenAddOns/PythonWrapper/ObjectPtr.h"

#if PY_MAJOR_VERSION < 3
#define PYTHON2
#endif

#define DEFINE_TYPE(type)           \
  template <>                       \
  bool is<type>(PyObject * obj);    \
  template <>                       \
  type get<type>(PyObject * obj);   \
  template <>                       \
  ObjectPtr set<type>(const type&); \
  static_assert(true, "")

namespace cepgen {
  namespace python {
    /// Import a Python module in a new reference-counted Python object
    ObjectPtr importModule(const std::string&);
    /// Define a Python module from a Python code in a new reference-counted Python object
    ObjectPtr defineModule(const std::string&, const std::string&);

    /// Check if a Python object holds a given C++ type
    template <typename T>
    bool is(PyObject* obj);
    /// Check if a Python object holds a given C++ type
    template <typename T>
    inline bool is(const ObjectPtr& obj) {
      return is<T>(obj.get());
    }

    /// Cast a Python object into a C++ type
    template <typename T>
    T get(PyObject* obj);
    /// Cast a Python object into a C++ type
    template <typename T>
    inline T get(const ObjectPtr& obj) {
      return get<T>(obj.get());
    }
    /// Build a new Python object from a C++ one
    template <typename T>
    ObjectPtr set(const T&);
    /// Build a new Python list from a C++ STL vector
    template <typename T>
    inline ObjectPtr set(const std::vector<T>& vec) {
      ObjectPtr list(PyList_New(vec.size()));
      for (size_t i = 0; i < vec.size(); ++i)
        PyList_SetItem(list.get(), i, set(vec.at(i)).release());
      return list;
    }

    /// Check if a Python object is compatible with a vector of uniform objects
    template <typename T>
    bool isVector(PyObject* obj);
    /// Check if a Python object is compatible with a vector of uniform objects
    template <typename T>
    inline bool isVector(const ObjectPtr& obj) {
      return isVector<T>(obj.get());
    }
    /// Retrieve a vector of objects, either from a Python list or tuple
    template <typename T>
    std::vector<T> getVector(PyObject* obj);
    /// Retrieve a vector of objects, either from a Python list or tuple
    template <typename T>
    inline std::vector<T> getVector(const ObjectPtr& obj) {
      return getVector<T>(obj.get());
    }

    /// Build a Python tuple from a C++ tuple
    template <typename... Args>
    ObjectPtr newTuple(const std::tuple<Args...>& c_tuple) {
      const auto tuple_size = sizeof...(Args);
      ObjectPtr tuple(PyTuple_New(tuple_size));
      Py_ssize_t i = 0;
      std::apply([&tuple, &i](
                     auto... vals) { ((PyTuple_SetItem(tuple.get(), i++, set<decltype(vals)>(vals).release())), ...); },
                 c_tuple);
      return tuple;
    }

    /// Build a Python tuple from a (uniform) vector of objects
    template <typename T>
    ObjectPtr newTuple(const std::vector<T>&);

    // type-specialised functions

    //--- long integers
    template <>
    bool is<long>(PyObject* obj);
    //--- unsigned long integers
    template <>
    unsigned long get<unsigned long>(PyObject* obj);
    //--- others
    DEFINE_TYPE(bool);
    DEFINE_TYPE(int);
    DEFINE_TYPE(double);
    DEFINE_TYPE(std::string);
    DEFINE_TYPE(Limits);
    DEFINE_TYPE(ParametersList);
  }  // namespace python
}  // namespace cepgen

#endif
