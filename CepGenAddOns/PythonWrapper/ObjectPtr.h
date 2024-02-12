/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2024  Laurent Forthomme
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

#ifndef CepGenAddOns_PythonWrapper_ObjectPtr_h
#define CepGenAddOns_PythonWrapper_ObjectPtr_h

#include <Python.h>

#include <memory>
#include <string>

namespace cepgen {
  namespace python {
    /// Common deleter for a PyObject
    struct ObjectPtrDeleter {
      void operator()(PyObject*);
    };
    /// Smart pointer to a Python object and its dereferencing operator
    struct ObjectPtr : public std::unique_ptr<PyObject, ObjectPtrDeleter> {
      using std::unique_ptr<PyObject, ObjectPtrDeleter>::unique_ptr;

      friend std::ostream& operator<<(std::ostream&, const ObjectPtr&);

      /// Retrieve the attribute from a python object
      ObjectPtr attribute(const std::string&) const;
    };
  }  // namespace python
}  // namespace cepgen

#endif
