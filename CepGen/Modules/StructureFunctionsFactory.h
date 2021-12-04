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

#ifndef CepGen_Modules_StructureFunctionsFactory_h
#define CepGen_Modules_StructureFunctionsFactory_h

#include "CepGen/Modules/ModuleFactory.h"

/// Add a structure functions definition to the list of handled parameterisation
#define REGISTER_STRFUN(id, name, obj)                                                             \
  namespace cepgen {                                                                               \
    struct BUILDERNM(name) {                                                                       \
      BUILDERNM(name)() { strfun::StructureFunctionsFactory::get().registerModule<obj>((int)id); } \
    };                                                                                             \
    static const BUILDERNM(name) gStrFun##name;                                                    \
  }
/// Add a structure functions definition (with its associated default parameters) to the list of handled parameterisation
#define REGISTER_STRFUN_PARAMS(id, name, obj, params)                                                      \
  namespace cepgen {                                                                                       \
    struct BUILDERNM(name) {                                                                               \
      BUILDERNM(name)() { strfun::StructureFunctionsFactory::get().registerModule<obj>((int)id, params); } \
    };                                                                                                     \
    static const BUILDERNM(name) gStrFun##name;                                                            \
  }

/// Add a sigma ratio definition to the list of handled parameterisation
#define SRBUILDERNM(id) Ratio##id##Builder
#define REGISTER_SIGRAT(id, name, obj)                                                        \
  namespace cepgen {                                                                          \
    struct SRBUILDERNM(name) {                                                                \
      SRBUILDERNM(name)() { sigrat::SigmaRatiosFactory::get().registerModule<obj>((int)id); } \
    };                                                                                        \
    static const SRBUILDERNM(name) gSigRat##name;                                             \
  }

/// Add a form factors definition to the list of handled parameterisation
#define REGISTER_FF_MODEL(name, obj)                                              \
  namespace cepgen {                                                              \
    namespace formfac {                                                           \
      struct BUILDERNM(obj) {                                                     \
        BUILDERNM(obj)() { FormFactorsFactory::get().registerModule<obj>(name); } \
      };                                                                          \
      static const BUILDERNM(obj) gFF##obj;                                       \
    }                                                                             \
  }

namespace cepgen {
  namespace strfun {
    class Parameterisation;
    /// A structure functions parameterisations factory
    DEFINE_FACTORY_INT(StructureFunctionsFactory,
                       Parameterisation,
                       "Nucleon structure functions parameterisations factory");
  }  // namespace strfun
  namespace sigrat {
    class Parameterisation;
    /// A sigma ratio parameterisations factory
    DEFINE_FACTORY_INT(SigmaRatiosFactory, Parameterisation, "Sigma L/T parameterisations factory");
  }  // namespace sigrat
  namespace formfac {
    class Parameterisation;
    /// A form factors parameterisations factory
    DEFINE_FACTORY_STR(FormFactorsFactory, Parameterisation, "Nucleon form factors factory");
    /// Standard dipole handler name
    static constexpr const char* gFFStandardDipoleHandler = "StandardDipole";
  }  // namespace formfac
}  // namespace cepgen

#endif
