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

#include <algorithm>
#include <string>

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/Filesystem.h"
#include "CepGen/Utils/String.h"
#include "CepGenAddOns/PythonWrapper/PythonUtils.h"

// clang-format off
#include <frameobject.h>
// clang-format on

#if PY_MAJOR_VERSION < 3
#define PYTHON2
#endif

namespace cepgen {
  namespace python {
    //------------------------------------------------------------------
    // Python API helpers
    //------------------------------------------------------------------

    std::string pythonPath(const std::string& file) {
      const auto dir = fs::path{file}.remove_filename();
      CG_DEBUG("PythonHandler") << "Adding {" << dir << "} to the default search paths.";
      utils::env::append("PYTHONPATH", dir);

      const auto filename = utils::replace_all(fs::path{file}.replace_extension("").string() /* remove the extension */,
                                               {{"../", ".."}, {"/", "."}});
      CG_DEBUG("PythonHandler") << "Python path: " << filename;
      return filename;
    }

    void error(const std::string& message) {
      throw CG_FATAL("PythonHandler:error").log([&message](auto& err) {
        PyObject *ptype = nullptr, *pvalue = nullptr, *ptraceback_obj = nullptr;
        // retrieve error indicator and clear it to handle ourself the error
        PyErr_Fetch(&ptype, &pvalue, &ptraceback_obj);
        PyErr_Clear();
        // ensure the objects retrieved are properly normalised and point to compatible objects
        PyErr_NormalizeException(&ptype, &pvalue, &ptraceback_obj);
        err << message;
        if (ptype != nullptr) {  // we can start the traceback
          err << "\nError: " << decode(PyObject_Str(pvalue));
          PyTracebackObject* ptraceback = (PyTracebackObject*)ptraceback_obj;
          std::string tabul = "↪ ";
          if (ptraceback != nullptr) {
            while (ptraceback->tb_next != nullptr) {
              PyFrameObject* pframe = ptraceback->tb_frame;
              if (pframe != nullptr) {
                int line = PyCode_Addr2Line(pframe->f_code, pframe->f_lasti);
                const auto filename = decode(pframe->f_code->co_filename), funcname = decode(pframe->f_code->co_name);
                err << utils::format(
                    "\n\t%s%s on %s (line %d)", tabul.c_str(), utils::boldify(funcname).c_str(), filename.c_str(), line);
              } else
                err << utils::format("\n\t%s issue in line %d", tabul.c_str(), ptraceback->tb_lineno);
              tabul = std::string("  ") + tabul;
              ptraceback = ptraceback->tb_next;
            }
          }
        }
        Py_Finalize();
      });
    }

    std::string decode(PyObject* obj) {
      if (!obj)
        return "(none)";
#ifdef PYTHON2
      return PyString_AsString(obj);
#else
      if (PyUnicode_Check(obj))
        return PyUnicode_AsUTF8(obj);
      if (PyBytes_Check(obj))
        return strdup(PyBytes_AS_STRING(obj));
      return "(none)";
#endif
    }

    ObjectPtr encode(const std::string& str) {
      ObjectPtr obj(PyUnicode_FromString(str.c_str()));  // new
      if (!obj)
        error("Failed to encode the following string:\n\t" + str);
      return obj;
    }

    PyObject* element(PyObject* obj, const std::string& key) {
      PyObject* pout = nullptr;
      {
        auto nink = encode(key);
        if (!nink)
          return pout;
        pout = PyDict_GetItem(obj, nink.get());  // borrowed
      }
      if (pout)
        CG_DEBUG("PythonHandler:element") << "retrieved " << pout->ob_type->tp_name << " element \"" << key << "\" "
                                          << "from " << obj->ob_type->tp_name << " object. "
                                          << "New reference count: " << pout->ob_refcnt;
      else
        CG_DEBUG("PythonHandler:element") << "did not retrieve a valid element \"" << key << "\"";
      return pout;
    }

    void fillParameter(PyObject* parent, const char* key, bool& out) {
      auto* pobj = element(parent, key);  // borrowed
      if (!pobj) {
        CG_DEBUG("PythonHandler") << "Failed to retrieve boolean object \"" << key << "\".";
        return;
      }
      try {
        out = (bool)get<int>(pobj);
      } catch (const LoggedException& e) {
        error(utils::format("Failed to retrieve boolean object \"%s\":\n\t", key) + e.message());
      }
    }

    void fillParameter(PyObject* parent, const char* key, int& out) {
      auto* pobj = element(parent, key);  // borrowed
      if (!pobj) {
        CG_DEBUG("PythonHandler") << "Failed to retrieve integer object \"" << key << "\".";
        return;
      }
      try {
        out = get<int>(pobj);
      } catch (const LoggedException& e) {
        error(utils::format("Failed to retrieve integer object \"%s\":\n\t", key) + e.message());
      }
    }

    void fillParameter(PyObject* parent, const char* key, unsigned long& out) {
      auto* pobj = element(parent, key);  // borrowed
      if (!pobj) {
        CG_DEBUG("PythonHandler") << "Failed to retrieve unsigned long integer object \"" << key << "\".";
        return;
      }
      try {
        out = get<unsigned long>(pobj);
      } catch (const LoggedException& e) {
        error(utils::format("Failed to retrieve unsigned long integer object \"%s\":\n\t", key) + e.message());
      }
    }

    void fillParameter(PyObject* parent, const char* key, unsigned int& out) {
      auto* pobj = element(parent, key);  // borrowed
      if (!pobj) {
        CG_DEBUG("PythonHandler") << "Failed to retrieve unsigned integer object \"" << key << "\".";
        return;
      }
      try {
        out = get<unsigned long>(pobj);
      } catch (const LoggedException& e) {
        error(utils::format("Failed to retrieve unsigned integer object \"%s\":\n\t", key) + e.message());
      }
    }

    void fillParameter(PyObject* parent, const char* key, double& out) {
      auto* pobj = element(parent, key);  // borrowed
      if (!pobj) {
        CG_DEBUG("PythonHandler") << "Failed to retrieve float object \"" << key << "\".";
        return;
      }
      try {
        out = get<double>(pobj);
      } catch (const LoggedException& e) {
        error(utils::format("Failed to retrieve float object \"%s\":\n\t", key) + e.message());
      }
    }

    void fillParameter(PyObject* parent, const char* key, std::string& out) {
      auto* pobj = element(parent, key);  // borrowed
      if (!pobj) {
        CG_DEBUG("PythonHandler") << "Failed to retrieve string object \"" << key << "\".";
        return;
      }
      try {
        out = get<std::string>(pobj);
      } catch (const LoggedException& e) {
        error(utils::format("Failed to retrieve string object \"%s\":\n\t", key) + e.message());
      }
    }

    void fillParameter(PyObject* obj, const char* key, Limits& out) {
      auto* pobj = element(obj, key);  // borrowed
      if (!pobj) {
        CG_DEBUG("PythonHandler") << "Failed to retrieve limits object \"" << key << "\".";
        return;
      }
      try {
        out = get<Limits>(pobj);
      } catch (const LoggedException& e) {
        error(utils::format("Failed to retrieve limits object \"%s\":\n\t", key) + e.message());
      }
    }

    void fillParameter(PyObject* parent, const char* key, std::vector<double>& out) {
      out.clear();
      auto* pobj = element(parent, key);  // borrowed
      if (!pobj) {
        CG_DEBUG("PythonHandler") << "Failed to retrieve floats collection object \"" << key << "\".";
        return;
      }
      try {
        out = getVector<double>(pobj);
      } catch (const LoggedException& e) {
        error(utils::format("Failed to retrieve floats collection object \"%s\":\n\t", key) + e.message());
      }
    }

    void fillParameter(PyObject* parent, const char* key, std::vector<std::string>& out) {
      out.clear();
      auto* pobj = element(parent, key);  // borrowed
      if (!pobj) {
        CG_DEBUG("PythonHandler") << "Failed to retrieve strings collection object \"" << key << "\".";
        return;
      }
      try {
        out = getVector<std::string>(pobj);
      } catch (const LoggedException& e) {
        error(utils::format("Failed to retrieve strings collection object \"%s\":\n\t", key) + e.message());
      }
    }

    void fillParameter(PyObject* parent, const char* key, std::vector<int>& out) {
      out.clear();
      auto* pobj = element(parent, key);  // borrowed
      if (!pobj) {
        CG_DEBUG("PythonHandler") << "Failed to retrieve integers collection object \"" << key << "\".";
        return;
      }
      try {
        out = getVector<int>(pobj);
      } catch (const LoggedException& e) {
        error(utils::format("Failed to retrieve integers collection object \"%s\":\n\t", key) + e.message());
      }
    }

    void fillParameter(PyObject* parent, const char* key, ParametersList& out) {
      auto* pobj = element(parent, key);  // borrowed
      if (!pobj) {
        CG_DEBUG("PythonHandler") << "Failed to retrieve parameters list object \"" << key << "\".";
        return;
      }
      try {
        out += get<ParametersList>(pobj);
      } catch (const LoggedException& e) {
        error(utils::format("Failed to retrieve parameters list object \"%s\":\n\t", key) + e.message());
      }
    }

    void fillParameter(PyObject* parent, const char* key, std::vector<ParametersList>& out) {
      out.clear();
      auto* pobj = element(parent, key);  // borrowed
      if (!pobj) {
        CG_DEBUG("PythonHandler") << "Failed to retrieve parameters list collection object \"" << key << "\".";
        return;
      }
      try {
        out = getVector<ParametersList>(pobj);
      } catch (const LoggedException& e) {
        error(utils::format("Failed to retrieve parameters list collection object \"%s\":\n\t", key) + e.message());
      }
    }
  }  // namespace python
}  // namespace cepgen
