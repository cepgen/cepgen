/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2021-2024  Laurent Forthomme
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

#ifndef CepGenAddOns_PythonWrapper_PythonConfigWriter_h
#define CepGenAddOns_PythonWrapper_PythonConfigWriter_h

#include <fstream>

namespace cepgen {
  class RunParameters;
  class ParametersDescription;
  namespace python {
    class PythonConfigWriter final {
    public:
      PythonConfigWriter(const std::string&);
      ~PythonConfigWriter();

      PythonConfigWriter& operator<<(const RunParameters&);
      PythonConfigWriter& operator<<(const ParametersDescription&);

    private:
      inline std::string offset(size_t num) { return std::string(num * tab_len_, ' '); }  ///< Compute a char-offset
      mutable std::ofstream file_;
      const size_t tab_len_{4};
    };
  }  // namespace python
}  // namespace cepgen

#endif
