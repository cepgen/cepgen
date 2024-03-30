/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2024  Laurent Forthomme
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

#ifndef CepGen_Modules_StructureFunctionsFactory_h
#define CepGen_Modules_StructureFunctionsFactory_h

#include "CepGen/Modules/ModuleFactory.h"

/// Add a structure functions definition to the list of handled parameterisation
#define REGISTER_STRFUN(name, obj)                                                                                    \
  namespace cepgen {                                                                                                  \
    namespace strfun {                                                                                                \
      struct BUILDERNM(obj) {                                                                                         \
        BUILDERNM(obj)() { StructureFunctionsFactory::get().addIndex(obj::index(), name).registerModule<obj>(name); } \
      };                                                                                                              \
      static const BUILDERNM(obj) gStrFun##obj;                                                                       \
    }                                                                                                                 \
  }                                                                                                                   \
  static_assert(true, "")

/// Add a sigma ratio definition to the list of handled parameterisation
#define REGISTER_SIGRAT(name, obj)                                                                             \
  namespace cepgen {                                                                                           \
    namespace sigrat {                                                                                         \
      struct BUILDERNM(obj) {                                                                                  \
        BUILDERNM(obj)() { SigmaRatiosFactory::get().addIndex(obj::index(), name).registerModule<obj>(name); } \
      };                                                                                                       \
      static const BUILDERNM(obj) gSigRat##obj;                                                                \
    }                                                                                                          \
  }                                                                                                            \
  static_assert(true, "")

/// Define a new factory instance for the definition of modules
#define DEFINE_HYBRID_FACTORY(name, obj_type, descr)                                               \
  DEFINE_FACTORY(Base##name, obj_type, descr);                                                     \
  class name : public Base##name {                                                                 \
  public:                                                                                          \
    using Base##name::Base##name;                                                                  \
    using Base##name::build;                                                                       \
    using Base##name::describeParameters;                                                          \
    static name& get();                                                                            \
    std::unique_ptr<obj_type> build(int, const ParametersList& = ParametersList()) const;          \
    ParametersDescription describeParameters(int, const ParametersList& = ParametersList()) const; \
    inline name& addIndex(int index, const std::string& mod_name) {                                \
      indices_[index] = mod_name;                                                                  \
      return *this;                                                                                \
    }                                                                                              \
                                                                                                   \
  private:                                                                                         \
    std::unordered_map<int, std::string> indices_;                                                 \
  };                                                                                               \
  static_assert(true, "")

namespace cepgen {
  namespace strfun {
    class Parameterisation;
  }
  /// A structure functions parameterisations factory
  DEFINE_HYBRID_FACTORY(StructureFunctionsFactory,
                        strfun::Parameterisation,
                        "Nucleon structure functions parameterisations factory");
  namespace sigrat {
    class Parameterisation;
  }
  /// A sigma ratio parameterisations factory
  DEFINE_HYBRID_FACTORY(SigmaRatiosFactory, sigrat::Parameterisation, "Sigma L/T parameterisations factory");
}  // namespace cepgen

#undef DEFINE_HYBRID_FACTORY

#endif
