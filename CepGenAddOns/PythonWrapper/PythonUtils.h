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

#ifndef CepGenAddOns_PythonWrapper_PythonUtils_h
#define CepGenAddOns_PythonWrapper_PythonUtils_h

#include "CepGenAddOns/PythonWrapper/PythonTypes.h"

namespace cepgen {
  namespace python {
    /// Translate a filename into a python-compatible path
    std::string pythonPath(const std::string&);
    /// Retrieve the element from a python dictionary
    PyObject* element(PyObject*, const std::string&);
    ObjectPtr getAttribute(PyObject*, const std::string&);
    inline ObjectPtr getAttribute(const ObjectPtr& obj, const std::string& attr) {
      return getAttribute(obj.get(), attr);
    }
    void setAttribute(const ObjectPtr& obj, const std::string& attr, const ObjectPtr& value);
    template <typename T>
    inline void setAttribute(const ObjectPtr& obj, const std::string& attr, const T& value) {
      setAttribute(obj, attr, set(value));
    }
    std::vector<std::wstring> info();
    /// Call a python function with a tuple of arguments
    ObjectPtr callArgs(const ObjectPtr& func, const ObjectPtr& args);
    /// Call the named method of a python object with an uncounted set of arguments
    template <typename... Args>
    inline ObjectPtr callMethod(const ObjectPtr& obj, const std::string& method, Args&&... args) {
      return call(getAttribute(obj, method), std::forward<Args>(args)...);
    }
    /// Call a python function with an uncounted set of arguments
    template <typename... Args>
    inline ObjectPtr call(const ObjectPtr& func, Args&&... args) {
      return callArgs(func, python::newTuple(std::make_tuple(std::forward<Args>(args)...)));  // new
    }

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
  }  // namespace python
}  // namespace cepgen

#endif
