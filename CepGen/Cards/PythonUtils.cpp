#include <algorithm>
#include <string>

#include "CepGen/Cards/PythonHandler.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Utils/String.h"

#include <frameobject.h>

#if PY_MAJOR_VERSION < 3
#define PYTHON2
#endif

namespace cepgen {
  namespace card {
    //------------------------------------------------------------------
    // Python API helpers
    //------------------------------------------------------------------

    std::string PythonHandler::pythonPath(const std::string& file) const {
      std::string s_filename = file;
      auto path = utils::split(s_filename, '/');
      if (path.size() > 1) {
        s_filename = *path.rbegin();
        path.pop_back();
        auto dir = utils::merge(path, "/");
        CG_DEBUG("PythonHandler") << "Adding \"" << dir << "\" to the default search paths.";
        setenv("PYTHONPATH", dir.c_str(), 1);
      }
      s_filename = s_filename.substr(0, s_filename.find_last_of("."));  // remove the extension
      utils::replace_all(s_filename, "../", "..");
      utils::replace_all(s_filename, "/", ".");
      CG_DEBUG("PythonHandler") << "Python path: " << s_filename;
      return s_filename;
    }

    void PythonHandler::throwPythonError(const std::string& message) const {
      PyObject *ptype = nullptr, *pvalue = nullptr, *ptraceback_obj = nullptr;
      // retrieve error indicator and clear it to handle ourself the error
      PyErr_Fetch(&ptype, &pvalue, &ptraceback_obj);
      PyErr_Clear();
      // ensure the objects retrieved are properly normalised and point to compatible objects
      PyErr_NormalizeException(&ptype, &pvalue, &ptraceback_obj);
      std::ostringstream oss;
      oss << message;
      if (ptype != nullptr) {  // we can start the traceback
        oss << "\n\tError: "
#ifdef PYTHON2
            << PyString_AsString(PyObject_Str(pvalue));  // deprecated in python v3+
#else
            << PyUnicode_AsUTF8(PyObject_Str(pvalue));
#endif
        PyTracebackObject* ptraceback = (PyTracebackObject*)ptraceback_obj;
        std::string tabul = "↪ ";
        if (ptraceback != nullptr) {
          while (ptraceback->tb_next != nullptr) {
            PyFrameObject* pframe = ptraceback->tb_frame;
            if (pframe != nullptr) {
              int line = PyCode_Addr2Line(pframe->f_code, pframe->f_lasti);
#ifdef PYTHON2
              const std::string filename = PyString_AsString(pframe->f_code->co_filename);
              const std::string = PyString_AsString(pframe->f_code->co_name);
#else
              const std::string filename = PyUnicode_AsUTF8(pframe->f_code->co_filename);
              const std::string funcname = PyUnicode_AsUTF8(pframe->f_code->co_name);
#endif
              oss << utils::format(
                  "\n\t%s%s on %s (line %d)", tabul.c_str(), utils::boldify(funcname).c_str(), filename, line);
            } else
              oss << utils::format("\n\t%s issue in line %d", tabul.c_str(), ptraceback->tb_lineno);
            tabul = std::string("  ") + tabul;
            ptraceback = ptraceback->tb_next;
          }
        }
      }
      Py_Finalize();
      throw CG_FATAL("PythonHandler:error") << oss.str();
    }

    PyObject* PythonHandler::encode(const char* str) const {
      PyObject* obj = PyUnicode_FromString(str);  // new
      if (!obj)
        throwPythonError(utils::format("Failed to encode the following string:\n\t%s", str));
      return obj;
    }

    PyObject* PythonHandler::element(PyObject* obj, const std::string& key) const {
      PyObject *pout = nullptr, *nink = encode(key.c_str());
      if (!nink)
        return pout;
      pout = PyDict_GetItem(obj, nink);  // borrowed
      Py_CLEAR(nink);
      if (pout)
        CG_DEBUG("PythonHandler:element") << "retrieved " << pout->ob_type->tp_name << " element \"" << key << "\" "
                                          << "from " << obj->ob_type->tp_name << " object\n\t"
                                          << "new reference count: " << pout->ob_refcnt;
      else
        CG_DEBUG("PythonHandler:element") << "did not retrieve a valid element \"" << key << "\"";
      return pout;
    }

    void PythonHandler::fillParameter(PyObject* parent, const char* key, bool& out) {
      PyObject* pobj = element(parent, key);  // borrowed
      if (!pobj) {
        CG_DEBUG("PythonHandler") << "Failed to retrieve boolean object \"" << key << "\".";
        return;
      }
      try {
        out = (bool)get<int>(pobj);
      } catch (const Exception& e) {
        throwPythonError(utils::format("Failed to retrieve boolean object \"%s\":\n\t%s", key, e.message().c_str()));
      }
    }

    void PythonHandler::fillParameter(PyObject* parent, const char* key, int& out) {
      PyObject* pobj = element(parent, key);  // borrowed
      if (!pobj) {
        CG_DEBUG("PythonHandler") << "Failed to retrieve integer object \"" << key << "\".";
        return;
      }
      try {
        out = get<int>(pobj);
      } catch (const Exception& e) {
        throwPythonError(utils::format("Failed to retrieve integer object \"%s\":\n\t%s", key, e.message().c_str()));
      }
    }

    void PythonHandler::fillParameter(PyObject* parent, const char* key, unsigned long& out) {
      PyObject* pobj = element(parent, key);  // borrowed
      if (!pobj) {
        CG_DEBUG("PythonHandler") << "Failed to retrieve unsigned long integer object \"" << key << "\".";
        return;
      }
      try {
        out = get<unsigned long>(pobj);
      } catch (const Exception& e) {
        throwPythonError(
            utils::format("Failed to retrieve unsigned long integer object \"%s\":\n\t%s", key, e.message().c_str()));
      }
    }

    void PythonHandler::fillParameter(PyObject* parent, const char* key, unsigned int& out) {
      PyObject* pobj = element(parent, key);  // borrowed
      if (!pobj) {
        CG_DEBUG("PythonHandler") << "Failed to retrieve unsigned integer object \"" << key << "\".";
        return;
      }
      try {
        out = get<unsigned long>(pobj);
      } catch (const Exception& e) {
        throwPythonError(
            utils::format("Failed to retrieve unsigned integer object \"%s\":\n\t%s", key, e.message().c_str()));
      }
    }

    void PythonHandler::fillParameter(PyObject* parent, const char* key, double& out) {
      PyObject* pobj = element(parent, key);  // borrowed
      if (!pobj) {
        CG_DEBUG("PythonHandler") << "Failed to retrieve float object \"" << key << "\".";
        return;
      }
      try {
        out = get<double>(pobj);
      } catch (const Exception& e) {
        throwPythonError(utils::format("Failed to retrieve float object \"%s\":\n\t%s", key, e.message().c_str()));
      }
    }

    void PythonHandler::fillParameter(PyObject* parent, const char* key, std::string& out) {
      PyObject* pobj = element(parent, key);  // borrowed
      if (!pobj) {
        CG_DEBUG("PythonHandler") << "Failed to retrieve string object \"" << key << "\".";
        return;
      }
      try {
        out = get<std::string>(pobj);
      } catch (const Exception& e) {
        throwPythonError(utils::format("Failed to retrieve string object \"%s\":\n\t%s", key, e.message().c_str()));
      }
    }

    void PythonHandler::fillParameter(PyObject* obj, const char* key, Limits& out) {
      PyObject* pobj = element(obj, key);  // borrowed
      if (!pobj) {
        CG_DEBUG("PythonHandler") << "Failed to retrieve limits object \"" << key << "\".";
        return;
      }
      try {
        out = get<Limits>(pobj);
      } catch (const Exception& e) {
        throwPythonError(utils::format("Failed to retrieve limits object \"%s\":\n\t%s", key, e.message().c_str()));
      }
    }

    void PythonHandler::fillParameter(PyObject* parent, const char* key, std::vector<double>& out) {
      out.clear();
      PyObject* pobj = element(parent, key);  // borrowed
      if (!pobj) {
        CG_DEBUG("PythonHandler") << "Failed to retrieve floats collection object \"" << key << "\".";
        return;
      }
      try {
        out = getVector<double>(pobj);
      } catch (const Exception& e) {
        throwPythonError(
            utils::format("Failed to retrieve floats collection object \"%s\":\n\t%s", key, e.message().c_str()));
      }
    }

    void PythonHandler::fillParameter(PyObject* parent, const char* key, std::vector<std::string>& out) {
      out.clear();
      PyObject* pobj = element(parent, key);  // borrowed
      if (!pobj) {
        CG_DEBUG("PythonHandler") << "Failed to retrieve strings collection object \"" << key << "\".";
        return;
      }
      try {
        out = getVector<std::string>(pobj);
      } catch (const Exception& e) {
        throwPythonError(
            utils::format("Failed to retrieve strings collection object \"%s\":\n\t%s", key, e.message().c_str()));
      }
    }

    void PythonHandler::fillParameter(PyObject* parent, const char* key, std::vector<int>& out) {
      out.clear();
      PyObject* pobj = element(parent, key);  // borrowed
      if (!pobj) {
        CG_DEBUG("PythonHandler") << "Failed to retrieve integers collection object \"" << key << "\".";
        return;
      }
      try {
        out = getVector<int>(pobj);
      } catch (const Exception& e) {
        throwPythonError(
            utils::format("Failed to retrieve integers collection object \"%s\":\n\t%s", key, e.message().c_str()));
      }
    }

    void PythonHandler::fillParameter(PyObject* parent, const char* key, ParametersList& out) {
      PyObject* pobj = element(parent, key);  // borrowed
      if (!pobj) {
        CG_DEBUG("PythonHandler") << "Failed to retrieve parameters list object \"" << key << "\".";
        return;
      }
      try {
        out += get<ParametersList>(pobj);
      } catch (const Exception& e) {
        throwPythonError(
            utils::format("Failed to retrieve parameters list object \"%s\":\n\t%s", key, e.message().c_str()));
      }
    }

    void PythonHandler::fillParameter(PyObject* parent, const char* key, std::vector<ParametersList>& out) {
      out.clear();
      PyObject* pobj = element(parent, key);  // borrowed
      if (!pobj) {
        CG_DEBUG("PythonHandler") << "Failed to retrieve parameters list collection object \"" << key << "\".";
        return;
      }
      try {
        out = getVector<ParametersList>(pobj);
      } catch (const Exception& e) {
        throwPythonError(utils::format(
            "Failed to retrieve parameters list collection object \"%s\":\n\t%s", key, e.message().c_str()));
      }
    }
  }  // namespace card
}  // namespace cepgen
