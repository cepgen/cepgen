/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2026  Laurent Forthomme
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

#include <fstream>

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/Filesystem.h"

namespace cepgen::utils {
  bool fileExists(const std::string& path) { return fs::exists(path); }

  std::string fileExtension(const std::string& file) { return fs::path(file).extension(); }

  std::string readFile(const std::string& filename) {
    if (auto input_file = std::ifstream(filename); input_file.good())
      return std::string(std::istreambuf_iterator(input_file), std::istreambuf_iterator<char>());
    throw CG_FATAL("readFile") << "Failed to open the file '" << filename << "' for reading.";
  }

  bool isWriteable(const std::string& path) {
    std::ofstream file(path);
    return file.good();
  }

  // Temporary working directory

  TemporaryDirectory::TemporaryDirectory(const std::string& base_dir, bool keep_files)
      : path_(std::filesystem::temp_directory_path() / base_dir), keep_files_(keep_files) {
    std::filesystem::create_directories(path_);
  }

  TemporaryDirectory::~TemporaryDirectory() {
    if (!keep_files_) {
      std::error_code ec;
      if (const auto num_files_removed = std::filesystem::remove_all(path_, ec);
          ec || num_files_removed == static_cast<std::uintmax_t>(-1))
        CG_ERROR("TemporaryDirectory") << "Failed to remove the temporary directory at " << path_ << ": "
                                       << ec.message() << ".";
      else
        CG_DEBUG("TemporaryDirectory") << "Successfully removed the temporary directory at " << path_ << ".";
    } else
      CG_DEBUG("TemporaryDirectory") << "Keeping the temporary directory structure at " << path_ << ".";
  }

  const fs::path& TemporaryDirectory::path() const { return path_; }
}  // namespace cepgen::utils
