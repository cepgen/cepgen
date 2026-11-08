/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2025  Laurent Forthomme
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

using namespace cepgen;

Steerable::Steerable(const ParametersList& params) { Steerable::setParameters(params); }

void Steerable::setParameters(const ParametersList& params) { params_ += params; }

std::string Steerable::steerPath(const std::string& key) const {
  if (const auto filename = steer<std::string>(key); !filename.empty()) {
    for (const auto& path : utils::env::searchPaths())
      if (const auto abs_path = fs::path(path) / filename; utils::fileExists(abs_path)) {
        CG_DEBUG("Steerable:steerPath") << "Found path for '" << key << "' at '" << abs_path << "'.";
        return abs_path;
      }
    return filename;
  }
  return {};
}

ParametersDescription Steerable::description() {
  auto desc = ParametersDescription("Steerable");
  desc.setDescription("Pure virtual base steerable object");
  return desc;
}
