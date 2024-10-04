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

#ifndef CepGenPython_Error_h
#define CepGenPython_Error_h

#include "CepGen/Core/Exception.h"
#include "CepGenPython/ObjectPtr.h"

#define PY_ERROR cepgen::python::Error(__FUNC__, __FILE__, __LINE__)

namespace cepgen::python {
  class Error final : public Exception {
  public:
    explicit Error(const char*, const char*, short) noexcept;

  private:
    PyObject* ptype_{nullptr};
    PyObject* pvalue_{nullptr};
    PyObject* ptraceback_obj_{nullptr};
  };
}  // namespace cepgen::python

#endif
