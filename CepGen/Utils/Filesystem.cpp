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

#include <unistd.h>

#include <fstream>

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/Filesystem.h"

namespace cepgen::utils {
  bool fileExists(const std::string& path) { return fs::exists(path); }

  std::string fileExtension(const std::string& file) { return fs::path(file).extension(); }

  std::string readFile(const std::string& filename) {
    if (auto input_file = std::ifstream(filename); input_file.good())
      return std::string(std::istreambuf_iterator<char>(input_file), std::istreambuf_iterator<char>());
    throw CG_FATAL("readFile") << "Failed to open the file '" << filename << "' for reading.";
  }

  bool isWriteable(const std::string& path) { return ::access(path.data(), W_OK) == 0; }
}  // namespace cepgen::utils
