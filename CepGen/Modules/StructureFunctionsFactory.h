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
#define REGISTER_SIGRAT(id, obj)                                                \
  namespace cepgen {                                                            \
    namespace sigrat {                                                          \
      struct BUILDERNM(obj) {                                                   \
        BUILDERNM(obj)() { SigmaRatiosFactory::get().registerModule<obj>(id); } \
      };                                                                        \
      static const BUILDERNM(obj) gSigRat##obj;                                 \
    }                                                                           \
  }                                                                             \
  static_assert(true, "")

namespace cepgen {
  namespace strfun {
    class Parameterisation;
  }
  /// A structure functions parameterisations base factory
  DEFINE_FACTORY(std::string,
                 BaseStructureFunctionsFactory,
                 strfun::Parameterisation,
                 "Nucleon structure functions parameterisations factory");
  /// A structure functions parameterisations factory
  class StructureFunctionsFactory : public BaseStructureFunctionsFactory {
  public:
    using BaseStructureFunctionsFactory::BaseStructureFunctionsFactory;
    using BaseStructureFunctionsFactory::build;
    using BaseStructureFunctionsFactory::describeParameters;
    static StructureFunctionsFactory& get();
    /// Build parameterisation from integer index
    std::unique_ptr<strfun::Parameterisation> build(int, const ParametersList& = ParametersList()) const;
    ParametersDescription describeParameters(int, const ParametersList& = ParametersList()) const;
    inline StructureFunctionsFactory& addIndex(int index, const std::string& name) {
      indices_[index] = name;
      return *this;
    }

  private:
    std::unordered_map<int, std::string> indices_;
  };
  namespace sigrat {
    class Parameterisation;
  }
  /// A sigma ratio parameterisations base factory
  DEFINE_FACTORY(std::string, BaseSigmaRatiosFactory, sigrat::Parameterisation, "Sigma L/T parameterisations factory");
  /// A sigma ratio parameterisations factory
  class SigmaRatiosFactory : public BaseSigmaRatiosFactory {
  public:
    using BaseSigmaRatiosFactory::BaseSigmaRatiosFactory;
    using BaseSigmaRatiosFactory::build;
    using BaseSigmaRatiosFactory::describeParameters;
    static SigmaRatiosFactory& get();
    /// Build parameterisation from integer index
    std::unique_ptr<sigrat::Parameterisation> build(int, const ParametersList& = ParametersList()) const;
    ParametersDescription describeParameters(int, const ParametersList& = ParametersList()) const;
    inline SigmaRatiosFactory& addIndex(int index, const std::string& name) {
      indices_[index] = name;
      return *this;
    }

  private:
    std::unordered_map<int, std::string> indices_;
  };
}  // namespace cepgen

#endif
