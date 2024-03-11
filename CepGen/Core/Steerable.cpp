/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2024  Laurent Forthomme
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

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/Steerable.h"
#include "CepGen/Utils/Environment.h"
#include "CepGen/Utils/Filesystem.h"

namespace cepgen {
  std::string Steerable::steerPath(const std::string& key) const {
    const auto fn = steer<std::string>(key);
    if (fn.empty())
      throw CG_ERROR("Steerable:steerPath") << "Trying to retrieve an empty path. Aborting.";
    for (const auto& path : utils::env::searchPaths())
      if (const auto abs_path = fs::path(path) / fn; utils::fileExists(abs_path)) {
        CG_DEBUG("Steerable:steerPath") << "Found path for '" << key << "' at '" << abs_path << "'.";
        return abs_path;
      }
    return fn;
  }
}  // namespace cepgen
