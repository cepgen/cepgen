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

#ifndef CepGen_Utils_Environment_h
#define CepGen_Utils_Environment_h

#include <string>
#include <vector>

namespace cepgen {
  namespace utils {
    /// All environment variable-related utilities
    namespace env {
      /// Retrieve a list of all search paths for external files
      std::vector<std::string> searchPaths();
      /// Get an environment variable
      std::string get(const std::string& var, const std::string& def = "");
      /// Set an environment variable
      void set(const std::string& var, const std::string& value);
      /// Add a value to an environment variable
      void append(const std::string& var, const std::string& value);
      /// Clear an environment variable
      void unset(const std::string& var);
    }  // namespace env
  }    // namespace utils
}  // namespace cepgen

#endif
