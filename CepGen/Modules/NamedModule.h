/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2021-2024  Laurent Forthomme
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

#include "CepGen/Core/SteeredObject.h"

namespace cepgen {
  /// Base runtime module object
  /// \tparam T Object type
  template <typename T>
  class NamedModule : public SteeredObject<T> {
  public:
    /// Build a module from its steering parameters
    explicit NamedModule(const ParametersList& params)
        : SteeredObject<T>(params), name_(SteeredObject<T>::steerName()) {}
    ~NamedModule() override = default;

    /// Describe all steering parameters for this module
    static ParametersDescription description() {
      auto desc = ParametersDescription();
      desc.setDescription("Named steerable module");
      return desc;
    }

    const std::string& name() const { return name_; }  ///< Module unique indexing name

  protected:
    const std::string name_;  ///< Module unique indexing name
  };
}  // namespace cepgen

#endif
