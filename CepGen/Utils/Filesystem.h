/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
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

#ifndef CepGen_Utils_Filesystem_h
#define CepGen_Utils_Filesystem_h
#if __cplusplus >= 201703L
#include <filesystem>
namespace fs = std::filesystem;
#elif __cplusplus >= 201103L
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#else
#error "*** no support for filesystem! ***"
#endif

#include <string>

namespace cepgen {
  namespace utils {
    /// Check if the file exists
    inline bool fileExists(const std::string& path) { return fs::exists(path); }
    /// Small utility to retrieve the extension of a filename
    inline std::string fileExtension(const std::string& file) { return fs::path(file).extension(); }
  }  //namespace utils
}  // namespace cepgen

#endif
