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

#ifndef CepGenPython_ObjectPtr_h
#define CepGenPython_ObjectPtr_h

#include <Python.h>

#if PY_MAJOR_VERSION < 3
#define PYTHON2
#endif

#include <memory>
#include <string>
#include <vector>

namespace cepgen::python {
  typedef std::unique_ptr<PyObject, void (*)(PyObject*)> PyObjectPtr;
  /// Smart pointer to a Python object and its dereferencing operator
  class ObjectPtr : public PyObjectPtr {
  public:
    explicit ObjectPtr(PyObject*, bool wrap_only = false);

    /// Wrap a PyObject without cleaning at the destructor
    static ObjectPtr wrap(PyObject*);
    /// Import a Python module in a new reference-counted Python object
    static ObjectPtr importModule(const std::string&);
    /// Define a Python module from a Python code in a new reference-counted Python object
    static ObjectPtr defineModule(const std::string&, const std::string&);

    /// Build a new Python object from a C++ one
    template <typename T>
    static ObjectPtr make(const T&);
    /// Build a new Python list from a C++ STL vector
    template <typename T>
    static inline ObjectPtr make(const std::vector<T>& vec) {
      ObjectPtr list(PyList_New(vec.size()));
      for (size_t i = 0; i < vec.size(); ++i)
        PyList_SetItem(list.get(), i, make(vec.at(i)).release());
      return list;
    }

    /// Check if a Python object holds a given C++ type
    template <typename T>
    bool is() const;
    /// Cast a Python object into a C++ type
    template <typename T>
    T value() const;
    /// Check if a Python object is compatible with a vector of uniform objects
    template <typename T>
    bool isVector() const;
    /// Retrieve a vector of objects, either from a Python list or tuple
    template <typename T>
    std::vector<T> vector() const;

    /// Build a Python tuple from a (uniform) vector of objects
    template <typename T>
    static ObjectPtr tupleFromVector(const std::vector<T>&);
    /// Build a Python tuple from a C++ tuple
    template <typename... Args>
    static inline ObjectPtr tuple(const std::tuple<Args...>& c_tuple) {
      const auto tuple_size = sizeof...(Args);
      ObjectPtr tuple(PyTuple_New(tuple_size));
      Py_ssize_t i = 0;
      std::apply(
          [&tuple, &i](auto... vals) {
            ((PyTuple_SetItem(tuple.get(), i++, make<decltype(vals)>(vals).release())), ...);
          },
          c_tuple);
      return tuple;
    }

    friend std::ostream& operator<<(std::ostream&, const ObjectPtr&);

    /// Call a python function with an uncounted set of arguments
    template <typename... Args>
    inline ObjectPtr operator()(Args&&... args) {
      return call(tuple(std::make_tuple(std::forward<Args>(args)...)));  // new
    }
    template <typename T>
    ObjectPtr call(const T&) const;                 ///< Call a python function with a single argument
    ObjectPtr call(const ObjectPtr&) const;         ///< Call a python function with a tuple of arguments
    ObjectPtr attribute(const std::string&) const;  ///< Retrieve the attribute from a python object
  };
}  // namespace cepgen::python

#endif
