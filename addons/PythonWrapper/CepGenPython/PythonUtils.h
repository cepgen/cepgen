/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2026  Laurent Forthomme
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

#ifndef CepGenPython_PythonUtils_h
#define CepGenPython_PythonUtils_h

#include <pybind11/pybind11.h>

namespace cepgen {
  class ParametersList;
}  // namespace cepgen

namespace cepgen::python {
  cepgen::ParametersList py_dict_to_plist(const pybind11::dict&);
  pybind11::dict plist_to_py_dict(const cepgen::ParametersList&);
}  // namespace cepgen::python

#endif
