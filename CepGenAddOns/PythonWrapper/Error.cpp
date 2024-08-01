/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2024  Laurent Forthomme
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
#include "CepGenAddOns/PythonWrapper/Utils.h"
#include <frameobject.h>
// clang-format on

#include "CepGen/Utils/String.h"

namespace cepgen::python {
  Error::Error(const char *origin, const char *file, short lineno) noexcept
      : Exception("", origin, Exception::Type::error, file, lineno) {
    // retrieve error indicator and clear it to handle ourself the error
    PyErr_Fetch(&ptype_, &pvalue_, &ptraceback_obj_);
    PyErr_Clear();
    // ensure the objects retrieved are properly normalised and point to compatible objects
    PyErr_NormalizeException(&ptype_, &pvalue_, &ptraceback_obj_);
    if (ptype_) {
      // we can start the traceback
      (*this) << "Error: " << ObjectPtr::wrap(PyObject_Str(pvalue_)).value<std::string>();
      if (auto mod_traceback = ObjectPtr::importModule("traceback"))
        if (auto fmt = mod_traceback.attribute("format_exception"); PyCallable_Check(fmt.get()))
          if (ObjectPtr pyth_val(PyObject_CallFunctionObjArgs(fmt.get(), ptype_, pvalue_, ptraceback_obj_, nullptr));
              pyth_val) {
            (*this) << "\n" << std::string(80, '.') << "\n";
            for (const auto &tb : pyth_val.vector<std::string>()) {
              size_t i = 0;
              std::string sep;
              for (const auto &err_line : utils::split(tb, '\n'))
                (*this) << sep << (i == 0 ? utils::boldify(err_line) : err_line), ++i, sep = "\n";
              (*this) << "\n";
            }
          }
      (*this) << std::string(80, '.') << "\n";
    }
  }
}  // namespace cepgen::python
