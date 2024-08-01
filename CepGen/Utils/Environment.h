/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2020-2024  Laurent Forthomme
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

/// All environment variable-related utilities
namespace cepgen::utils::env {
  std::vector<std::string> searchPaths();  ///< Retrieve a list of all search paths for external files
  std::string get(const std::string& var, const std::string& def = "");  ///< Get an environment variable
  void set(const std::string& var, const std::string& value);            ///< Set an environment variable
  void append(const std::string& var, const std::string& value);         ///< Add a value to an environment variable
  void unset(const std::string& var);                                    ///< Clear an environment variable
}  // namespace cepgen::utils::env

#endif
