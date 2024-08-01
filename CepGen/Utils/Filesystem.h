/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2024  Laurent Forthomme
 *                2023       Dmitri Konstantinov
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

#if defined(__has_include)
#if __has_include(<filesystem>)
#include <filesystem>
namespace fs = std::filesystem;
#elif __has_include(<experimental/filesystem>)
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#else
#error "*** no support for filesystem! ***"
#endif
#else
#error "*** no support for __has_include! ***"
#endif

#include <string>

namespace cepgen {
  namespace utils {
    bool fileExists(const std::string&);                 ///< Check if the file exists
    std::string fileExtension(const std::string& file);  ///< Small utility to retrieve the extension of a filename
    std::string readFile(const std::string&);            ///< Read the content of a file into a string buffer
    bool isWriteable(const std::string&);                ///< Check if path can be accessed for writing
  }  //namespace utils
}  // namespace cepgen

#endif
