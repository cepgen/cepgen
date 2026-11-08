/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2025  Laurent Forthomme
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

/// Add a structure functions definition to the list of handled parameterisations
#define REGISTER_STRFUN(name, id, obj)                                                                       \
  namespace cepgen::strfun {                                                                                 \
    struct BUILDER_NAME(obj) {                                                                               \
      BUILDER_NAME(obj)() { StructureFunctionsFactory::get().addIndex(id, name).registerModule<obj>(name); } \
    };                                                                                                       \
    static const BUILDER_NAME(obj) gStrFun##obj;                                                             \
  }                                                                                                          \
  static_assert(true, "")

/// Add a longitudinal/transverse cross-sections ratio definition to the list of handled parameterisations
#define REGISTER_SIGMA_RATIO(name, id, obj)                                                           \
  namespace cepgen::sigrat {                                                                          \
    struct BUILDER_NAME(obj) {                                                                        \
      BUILDER_NAME(obj)() { SigmaRatiosFactory::get().addIndex(id, name).registerModule<obj>(name); } \
    };                                                                                                \
    static const BUILDER_NAME(obj) gSigRat##obj;                                                      \
  }                                                                                                   \
  static_assert(true, "")

namespace cepgen::strfun {
  class Parameterisation;
}  // namespace cepgen::strfun
namespace cepgen::sigrat {
  class Parameterisation;
}  // namespace cepgen::sigrat

namespace cepgen {
  /// A structure functions parameterisations factory
  DEFINE_FACTORY(StructureFunctionsFactory,
                 strfun::Parameterisation,
                 "Nucleon structure functions parameterisations factory");
  /// A sigma ratio parameterisations factory
  DEFINE_FACTORY(SigmaRatiosFactory, sigrat::Parameterisation, "Sigma L/T parameterisations factory");
}  // namespace cepgen

#endif
