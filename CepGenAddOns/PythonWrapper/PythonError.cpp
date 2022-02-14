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

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/String.h"
#include "CepGenAddOns/PythonWrapper/PythonError.h"
#include "CepGenAddOns/PythonWrapper/PythonUtils.h"

// clang-format off
#include <frameobject.h>
// clang-format on

namespace cepgen {
  namespace python {
    Error::Error(const std::string& origin, const std::string& file, short lineno)
        : Exception(origin.c_str(), "Python::error", Exception::Type::error, file.c_str(), lineno) {
      // retrieve error indicator and clear it to handle ourself the error
      PyErr_Fetch(&ptype_, &pvalue_, &ptraceback_obj_);
      PyErr_Clear();
      // ensure the objects retrieved are properly normalised and point to compatible objects
      PyErr_NormalizeException(&ptype_, &pvalue_, &ptraceback_obj_);
    }

    Error::~Error() {
      if (ptype_) {  // we can start the traceback
        (*this) << "\nError: " << decode(PyObject_Str(pvalue_));
        auto* ptraceback = (PyTracebackObject*)ptraceback_obj_;
        std::string tabul = "↪ ";
        if (ptraceback) {
          while (ptraceback->tb_next) {
            auto* pframe = ptraceback->tb_frame;
            if (pframe) {
              int line = PyCode_Addr2Line(pframe->f_code, pframe->f_lasti);
              const auto filename = decode(pframe->f_code->co_filename), funcname = decode(pframe->f_code->co_name);
              (*this) << utils::format(
                  "\n\t%s%s on %s (line %d)", tabul.c_str(), utils::boldify(funcname).c_str(), filename.c_str(), line);
            } else
              (*this) << utils::format("\n\t%s issue in line %d", tabul.c_str(), ptraceback->tb_lineno);
            tabul = std::string("  ") + tabul;
            ptraceback = ptraceback->tb_next;
          }
        }
      }
      Py_Finalize();
      CG_LOG << "hahah";
    }
  }  // namespace python
}  // namespace cepgen
