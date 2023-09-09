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
#include "CepGenAddOns/PythonWrapper/Error.h"
#include "CepGenAddOns/PythonWrapper/PythonUtils.h"
#include <frameobject.h>
// clang-format on

#if PY_VERSION_HEX < 0x030b0000
#include "CepGen/Utils/String.h"
#endif

namespace cepgen {
  namespace python {
    Error::Error(const char* origin, const char* file, short lineno) noexcept
        : Exception("", origin, Exception::Type::error, file, lineno) {
      // retrieve error indicator and clear it to handle ourself the error
      PyErr_Fetch(&ptype_, &pvalue_, &ptraceback_obj_);
      PyErr_Clear();
      // ensure the objects retrieved are properly normalised and point to compatible objects
      PyErr_NormalizeException(&ptype_, &pvalue_, &ptraceback_obj_);
      if (ptype_) {
        // we can start the traceback
        (*this) << "Error: " << get<std::string>(PyObject_Str(pvalue_));
        if (auto* ptraceback = reinterpret_cast<PyTracebackObject*>(ptraceback_obj_)) {
          const std::string arr = "↪ ";
          std::string tabul;
          while (ptraceback->tb_next) {
            (*this) << "\n\t" << tabul << arr;
#if PY_VERSION_HEX < 0x030b0000
            if (PyFrameObject* pframe = ptraceback->tb_frame)
              (*this) << utils::boldify(get<std::string>(pframe->f_code->co_name)) << " on "
                      << get<std::string>(pframe->f_code->co_filename) << " (line "
                      << PyCode_Addr2Line(pframe->f_code, pframe->f_lasti) << ")";
            else
#endif
              (*this) << " issue on line " << ptraceback->tb_lineno;
            tabul += "  ";
            ptraceback = ptraceback->tb_next;
          }
        }
        (*this) << "\n" << std::string(80, '.') << "\n";
      }
    }
  }  // namespace python
}  // namespace cepgen
