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

// clang-format off
#include "CepGenAddOns/PythonWrapper/PythonError.h"
#include "CepGenAddOns/PythonWrapper/PythonUtils.h"
// clang-format on

#include <algorithm>

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/Filesystem.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  namespace python {
    //------------------------------------------------------------------
    // Python API helpers
    //------------------------------------------------------------------

    Environment::Environment() {
      Py_InitializeEx(1);
      if (!initialised())
        throw CG_FATAL("Python:Environment") << "Failed to initialise the Python environment!";
    }

    Environment::~Environment() {
      if (!initialised())
        CG_FATAL("Python:Environment")
            << "Python environment is set to be finalised while it was not initialised in the first place.";
      Py_Finalize();
    }

    bool Environment::initialised() { return Py_IsInitialized(); }

    std::string pythonPath(const std::string& file) {
      const auto dir = fs::path{file}.remove_filename();
      CG_DEBUG("Python") << "Adding {" << dir << "} to the default search paths.";
      utils::env::append("PYTHONPATH", dir);

      const auto filename = utils::replace_all(fs::path{file}.replace_extension("").string() /* remove the extension */,
                                               {{"../", ".."}, {"/", "."}});
      CG_DEBUG("Python") << "Python path: " << filename;
      return filename;
    }

    PyObject* element(PyObject* obj, const std::string& key) {
      auto* pout = PyDict_GetItem(obj, set(key).release());  // borrowed
      if (pout)
        CG_DEBUG("Python:element") << "retrieved " << pout->ob_type->tp_name << " element \"" << key << "\" "
                                   << "from " << obj->ob_type->tp_name << " object. "
                                   << "New reference count: " << pout->ob_refcnt;
      else
        CG_DEBUG("Python:element") << "did not retrieve a valid element \"" << key << "\"";
      return pout;
    }

    ObjectPtr getAttribute(PyObject* obj, const std::string& attr) {
      if (PyObject_HasAttrString(obj, attr.c_str()) != 1)
        return ObjectPtr(nullptr);
      return ObjectPtr(PyObject_GetAttrString(obj, attr.c_str()));  // new
    }

    void fillParameter(PyObject* parent, const char* key, bool& out) {
      auto* pobj = element(parent, key);  // borrowed
      if (!pobj) {
        CG_DEBUG("Python") << "Failed to retrieve boolean object \"" << key << "\".";
        return;
      }
      try {
        out = (bool)get<int>(pobj);
      } catch (const Exception& e) {
        PY_ERROR << "Failed to retrieve boolean object \"" << key << "\":\n\t" << e.message();
      }
    }

    void fillParameter(PyObject* parent, const char* key, int& out) {
      auto* pobj = element(parent, key);  // borrowed
      if (!pobj) {
        CG_DEBUG("Python") << "Failed to retrieve integer object \"" << key << "\".";
        return;
      }
      try {
        out = get<int>(pobj);
      } catch (const Exception& e) {
        PY_ERROR << "Failed to retrieve integer object \"" << key << "\":\n\t" << e.message();
      }
    }

    void fillParameter(PyObject* parent, const char* key, unsigned long& out) {
      auto* pobj = element(parent, key);  // borrowed
      if (!pobj) {
        CG_DEBUG("Python") << "Failed to retrieve unsigned long integer object \"" << key << "\".";
        return;
      }
      try {
        out = get<unsigned long>(pobj);
      } catch (const Exception& e) {
        PY_ERROR << "Failed to retrieve unsigned long integer object \"" << key << "\":\n\t" << e.message();
      }
    }

    void fillParameter(PyObject* parent, const char* key, unsigned int& out) {
      auto* pobj = element(parent, key);  // borrowed
      if (!pobj) {
        CG_DEBUG("Python") << "Failed to retrieve unsigned integer object \"" << key << "\".";
        return;
      }
      try {
        out = get<unsigned long>(pobj);
      } catch (const Exception& e) {
        PY_ERROR << "Failed to retrieve unsigned integer object \"" << key << "\":\n\t" << e.message();
      }
    }

    void fillParameter(PyObject* parent, const char* key, double& out) {
      auto* pobj = element(parent, key);  // borrowed
      if (!pobj) {
        CG_DEBUG("Python") << "Failed to retrieve float object \"" << key << "\".";
        return;
      }
      try {
        out = get<double>(pobj);
      } catch (const Exception& e) {
        PY_ERROR << "Failed to retrieve float object \"" << key << "\":\n\t" << e.message();
      }
    }

    void fillParameter(PyObject* parent, const char* key, std::string& out) {
      auto* pobj = element(parent, key);  // borrowed
      if (!pobj) {
        CG_DEBUG("Python") << "Failed to retrieve string object \"" << key << "\".";
        return;
      }
      try {
        out = get<std::string>(pobj);
      } catch (const Exception& e) {
        PY_ERROR << "Failed to retrieve string object \"" << key << "\":\n\t" << e.message();
      }
    }

    void fillParameter(PyObject* obj, const char* key, Limits& out) {
      auto* pobj = element(obj, key);  // borrowed
      if (!pobj) {
        CG_DEBUG("Python") << "Failed to retrieve limits object \"" << key << "\".";
        return;
      }
      try {
        out = get<Limits>(pobj);
      } catch (const Exception& e) {
        PY_ERROR << "Failed to retrieve limits object \"" << key << "\":\n\t" << e.message();
      }
    }

    void fillParameter(PyObject* parent, const char* key, std::vector<double>& out) {
      out.clear();
      auto* pobj = element(parent, key);  // borrowed
      if (!pobj) {
        CG_DEBUG("Python") << "Failed to retrieve floats collection object \"" << key << "\".";
        return;
      }
      try {
        out = getVector<double>(pobj);
      } catch (const Exception& e) {
        PY_ERROR << "Failed to retrieve floats collection object \"" << key << "\":\n\t" << e.message();
      }
    }

    void fillParameter(PyObject* parent, const char* key, std::vector<std::string>& out) {
      out.clear();
      auto* pobj = element(parent, key);  // borrowed
      if (!pobj) {
        CG_DEBUG("Python") << "Failed to retrieve strings collection object \"" << key << "\".";
        return;
      }
      try {
        out = getVector<std::string>(pobj);
      } catch (const Exception& e) {
        PY_ERROR << "Failed to retrieve strings collection object \"" << key << "\":\n\t" << e.message();
      }
    }

    void fillParameter(PyObject* parent, const char* key, std::vector<int>& out) {
      out.clear();
      auto* pobj = element(parent, key);  // borrowed
      if (!pobj) {
        CG_DEBUG("Python") << "Failed to retrieve integers collection object \"" << key << "\".";
        return;
      }
      try {
        out = getVector<int>(pobj);
      } catch (const Exception& e) {
        PY_ERROR << "Failed to retrieve integers collection object \"" << key << "\":\n\t" << e.message();
      }
    }

    void fillParameter(PyObject* parent, const char* key, ParametersList& out) {
      auto* pobj = element(parent, key);  // borrowed
      if (!pobj) {
        CG_DEBUG("Python") << "Failed to retrieve parameters list object \"" << key << "\".";
        return;
      }
      try {
        out += get<ParametersList>(pobj);
      } catch (const Exception& e) {
        PY_ERROR << "Failed to retrieve parameters list object \"" << key << "\":\n\t" << e.message();
      }
    }

    void fillParameter(PyObject* parent, const char* key, std::vector<ParametersList>& out) {
      out.clear();
      auto* pobj = element(parent, key);  // borrowed
      if (!pobj) {
        CG_DEBUG("Python") << "Failed to retrieve parameters list collection object \"" << key << "\".";
        return;
      }
      try {
        out = getVector<ParametersList>(pobj);
      } catch (const Exception& e) {
        PY_ERROR << "Failed to retrieve parameters list collection object \"" << key << "\":\n\t" << e.message();
      }
    }
  }  // namespace python
}  // namespace cepgen
