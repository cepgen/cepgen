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

#include "CepGen/Core/Steerable.h"

namespace cepgen {
  /// Base runtime module object
  /// \tparam I Indexing type for runtime factory builder registration
  template <typename I = std::string>
  class NamedModule : public Steerable {
  public:
    /// Build a module from its steering parameters
    explicit NamedModule(const ParametersList& params) : Steerable(params), name_(params.name<I>()) {}
    virtual ~NamedModule() = default;

    /// Describe all steering parameters for this module
    static ParametersDescription description() {
      auto desc = ParametersDescription();
      desc.setDescription("Named steerable module");
      return desc;
    }

    /// Module unique indexing name
    const I& name() const { return name_; }

  protected:
    /// Module unique indeing name
    const I name_;
  };
}  // namespace cepgen

#endif
