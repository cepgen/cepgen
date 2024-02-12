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

#include "CepGen/Utils/Message.h"
#include "CepGenAddOns/PythonWrapper/ObjectPtr.h"
#include "CepGenAddOns/PythonWrapper/PythonTypes.h"

namespace cepgen {
  namespace python {
    void ObjectPtrDeleter::operator()(PyObject* obj) {
      CG_DEBUG("Python:ObjectPtrDeleter").log([&obj](auto& log) {
        log << "Destroying object at addr 0x" << obj << " (";
#if PY_VERSION_HEX >= 0x03110000
        if (auto* type = Py_TYPE(obj); type)
          log << "type: " << get<std::string>(PyType_GetName(type)) << ", ";
#endif
        log << "reference count: " << Py_REFCNT(obj) << ")";
      });
      Py_DECREF(obj);
    }

    std::ostream& operator<<(std::ostream& os, const ObjectPtr& ptr) {
      os << "PyObject{";
      if (auto repr = ObjectPtr(PyObject_Str(ptr.get())); repr)  // new
        os << get<std::string>(repr);
      return os << "}";
    }

    ObjectPtr ObjectPtr::attribute(const std::string& attr) const {
      if (PyObject_HasAttrString(get(), attr.c_str()) != 1)
        return ObjectPtr(nullptr);
      return ObjectPtr(PyObject_GetAttrString(get(), attr.c_str()));  // new
    }
  }  // namespace python
}  // namespace cepgen
