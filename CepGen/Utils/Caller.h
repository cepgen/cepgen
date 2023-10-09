/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023  Laurent Forthomme
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

#ifndef CepGen_Utils_Caller_h
#define CepGen_Utils_Caller_h

#include <string>
#include <vector>

namespace cepgen {
  namespace utils {
    /// External command piping utility
    class Caller {
    public:
      Caller() {}
      /// Start a logged call command
      /// \param[in] commands Command path for the session
      static std::string call(const std::vector<std::string>& commands);
      /// Start a logged call command
      /// \param[in] command Command path for the session
      static std::string call(const std::string& command);
    };
  }  // namespace utils
}  // namespace cepgen

#endif
