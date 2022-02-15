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

#ifndef CepGenAddOns_PythonWrapper_PythonError_h
#define CepGenAddOns_PythonWrapper_PythonError_h

#include <Python.h>

#include <string>

#define PY_ERROR cepgen::python::Error(__FUNC__, __FILE__, __LINE__)

namespace cepgen {
  namespace python {
    class Error final : public Exception {
    public:
      explicit Error(const std::string&, const std::string&, short);
      ~Error();

    private:
      PyObject* ptype_{nullptr};
      PyObject* pvalue_{nullptr};
      PyObject* ptraceback_obj_{nullptr};
    };
  }  // namespace python
}  // namespace cepgen

#endif
