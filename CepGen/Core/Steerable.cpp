/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2023  Laurent Forthomme
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

#include "CepGen/Core/Steerable.h"
#include "CepGen/Utils/Environment.h"
#include "CepGen/Utils/Filesystem.h"
#include "CepGen/Utils/Message.h"

namespace cepgen {
  std::string Steerable::steerPath(const std::string& key) const {
    const auto& fn = steer<std::string>(key);
    if (fn.empty())
      return fn;
    for (const auto& path : utils::env::searchPaths())
      if (utils::fileExists(fs::path(path) / fn)) {
        CG_DEBUG("Steerable:steerPath") << "Found path for '" << key << "' at '" << fs::path(path) / fn << "'.";
        return fs::path(path) / fn;
      }
    return fn;
  }
}  // namespace cepgen
