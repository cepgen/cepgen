/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023  Laurent Forthomme
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

#ifndef CepGen_Modules_PartonFluxFactory_h
#define CepGen_Modules_PartonFluxFactory_h

#include "CepGen/Modules/ModuleFactory.h"

/** \file */

/// Add a generic parton flux evaluator builder definition
#define REGISTER_FLUX(name, obj)                                               \
  namespace cepgen {                                                           \
    struct BUILDERNM(obj) {                                                    \
      BUILDERNM(obj)() { PartonFluxFactory::get().registerModule<obj>(name); } \
    };                                                                         \
    static const BUILDERNM(obj) gPartonFlux##obj;                              \
  }                                                                            \
  static_assert(true, "")

/// Add a collinear flux definition to the list of handled parameterisation
#define REGISTER_COLLFLUX(name, obj)                                              \
  namespace cepgen {                                                              \
    struct BUILDERNM(obj) {                                                       \
      BUILDERNM(obj)() { CollinearFluxFactory::get().registerModule<obj>(name); } \
    };                                                                            \
    static const BUILDERNM(obj) gCollFlux##obj;                                   \
  }

namespace cepgen {
  class PartonFlux;
  namespace collflux {
    class Parameterisation;
  }
  /// A generic flux parameterisations factory
  DEFINE_FACTORY_STR(PartonFluxFactory, PartonFlux, "Parton flux estimators factory");
  /// A collinear flux parameterisations factory
  DEFINE_FACTORY_STR(CollinearFluxFactory, collflux::Parameterisation, "Collinear parton flux estimators factory");
}  // namespace cepgen

#endif
