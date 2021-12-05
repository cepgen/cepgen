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

#ifndef CepGen_Modules_NamedModule_h
#define CepGen_Modules_NamedModule_h

#include "CepGen/Utils/ParametersDescription.h"

namespace cepgen {
  /// Base runtime module object
  template <typename T = std::string>
  class NamedModule {
  public:
    /// Build a module from its steering parameters
    explicit NamedModule(const ParametersList& params) : params_(params), name_(params.name<T>()) {}
    virtual ~NamedModule() = default;

    /// Module unique name
    const T& name() const { return name_; }
    /// Description of all module parameters
    static inline ParametersDescription description() {
      auto desc = ParametersDescription();
      desc.setDescription("Unnamed module");
      return desc;
    }
    /// Module user-defined parameters
    inline const ParametersList& parameters() const { return params_; }

  protected:
    /// Set of parameters to steer this output module
    const ParametersList params_;
    /// Module unique name
    const T name_;
  };
}  // namespace cepgen

#endif
