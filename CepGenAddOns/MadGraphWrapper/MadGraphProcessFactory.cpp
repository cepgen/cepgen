/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2024  Laurent Forthomme
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
#include "CepGen/Modules/ModuleFactory.h"
#include "CepGenAddOns/MadGraphWrapper/MadGraphProcess.h"

namespace cepgen {
  template class ModuleFactory<MadGraphProcess, std::string>;

  template <>
  ModuleFactory<MadGraphProcess, std::string>::ModuleFactory(const std::string& desc) : description_(desc) {}

  template <>
  std::vector<std::string> ModuleFactory<MadGraphProcess, std::string>::modules() const {
    std::vector<std::string> out;
    std::transform(map_.begin(), map_.end(), std::back_inserter(out), [](const auto& val) { return val.first; });
    std::sort(out.begin(), out.end());
    return out;
  }

  template <>
  std::unique_ptr<MadGraphProcess> ModuleFactory<MadGraphProcess, std::string>::build(
      const std::string& mod_name, const ParametersList& params) const {
    if (map_.count(mod_name) == 0)
      throw CG_FATAL("ModuleFactory") << description_ << " failed to build a MadGraph5 process with name '" << mod_name
                                      << "'. Registered modules: " << modules() << ".";
    return map_.at(mod_name)(params_map_.at(mod_name).validate(params));
  }
}  // namespace cepgen
