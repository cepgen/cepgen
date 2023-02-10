/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2023  Laurent Forthomme
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
#define REGISTER_STRFUN(id, name, obj)                                                     \
  namespace cepgen {                                                                       \
    struct BUILDERNM(name) {                                                               \
      BUILDERNM(name)() { StructureFunctionsFactory::get().registerModule<obj>((int)id); } \
    };                                                                                     \
    static const BUILDERNM(name) gStrFun##name;                                            \
  }                                                                                        \
  static_assert(true, "")

/// Add a sigma ratio definition to the list of handled parameterisation
#define SRBUILDERNM(id) Ratio##id##Builder
#define REGISTER_SIGRAT(id, name, obj)                                                \
  namespace cepgen {                                                                  \
    struct SRBUILDERNM(name) {                                                        \
      SRBUILDERNM(name)() { SigmaRatiosFactory::get().registerModule<obj>((int)id); } \
    };                                                                                \
    static const SRBUILDERNM(name) gSigRat##name;                                     \
  }                                                                                   \
  static_assert(true, "")

namespace cepgen {
  namespace strfun {
    class Parameterisation;
  }
  /// A structure functions parameterisations factory
  DEFINE_FACTORY_INT(StructureFunctionsFactory,
                     strfun::Parameterisation,
                     "Nucleon structure functions parameterisations factory");
  namespace sigrat {
    class Parameterisation;
  }
  /// A sigma ratio parameterisations factory
  DEFINE_FACTORY_INT(SigmaRatiosFactory, sigrat::Parameterisation, "Sigma L/T parameterisations factory");
}  // namespace cepgen

#endif
