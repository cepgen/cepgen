/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024  Laurent Forthomme
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
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

#define IMPLEMENT_HYBRID_FACTORY(name, obj_type)                                                                     \
  name& name::get() {                                                                                                \
    static name instance;                                                                                            \
    return instance;                                                                                                 \
  }                                                                                                                  \
  std::unique_ptr<obj_type> name::build(int index, const ParametersList& params) const {                             \
    if (indices_.count(index) > 0)                                                                                   \
      return name::build(indices_.at(index), params);                                                                \
    const auto& mod_names = modules();                                                                               \
    if (const auto str_index = std::to_string(index);                                                                \
        std::find(mod_names.begin(), mod_names.end(), str_index) != mod_names.end())                                 \
      return name::build(str_index, params);                                                                         \
    throw CG_FATAL("HybridFactory") << description() << " failed to build a module with index '" << index << "'. \n" \
                                    << "Registered indices: " << indices_ << ".";                                    \
  }                                                                                                                  \
  ParametersDescription name::describeParameters(int index, const ParametersList& params) const {                    \
    if (indices_.count(index) > 0)                                                                                   \
      return name::describeParameters(indices_.at(index), params);                                                   \
    const auto& mod_names = modules();                                                                               \
    if (const auto str_index = std::to_string(index);                                                                \
        std::find(mod_names.begin(), mod_names.end(), str_index) != mod_names.end())                                 \
      return name::describeParameters(str_index, params);                                                            \
    throw CG_FATAL("HybridFactory") << "No parameters description were found for module index '" << index << "'.\n"  \
                                    << "Registered modules: " << indices_ << ".";                                    \
  }                                                                                                                  \
  static_assert(true, "")

namespace cepgen {
  IMPLEMENT_HYBRID_FACTORY(StructureFunctionsFactory, strfun::Parameterisation);
  IMPLEMENT_HYBRID_FACTORY(SigmaRatiosFactory, sigrat::Parameterisation);
}  // namespace cepgen
