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

#ifndef CepGen_Modules_CouplingFactory_h
#define CepGen_Modules_CouplingFactory_h

#include "CepGen/Modules/ModuleFactory.h"

/** \file */

/// Add an electromagnetic coupling evolution algorithm
#define REGISTER_ALPHAEM_MODULE(name, obj)                                  \
  namespace cepgen {                                                        \
    struct BUILDERNM(obj) {                                                 \
      BUILDERNM(obj)() { AlphaEMFactory::get().registerModule<obj>(name); } \
    };                                                                      \
    static const BUILDERNM(obj) gAlphaEM##obj;                              \
  }                                                                         \
  static_assert(true, "")

/// Add a strong coupling evolution algorithm
#define REGISTER_ALPHAS_MODULE(name, obj)                                  \
  namespace cepgen {                                                       \
    struct BUILDERNM(obj) {                                                \
      BUILDERNM(obj)() { AlphaSFactory::get().registerModule<obj>(name); } \
    };                                                                     \
    static const BUILDERNM(obj) gAlphaS##obj;                              \
  }                                                                        \
  static_assert(true, "")

namespace cepgen {
  class Coupling;
  /// An electromagnetic coupling evolution algorithms factory
  DEFINE_FACTORY_STR(AlphaEMFactory, Coupling, "Electromagnetic coupling evolution factory");
  /// A strong coupling evolution algorithms factory
  DEFINE_FACTORY_STR(AlphaSFactory, Coupling, "Strong coupling evolution factory");
}  // namespace cepgen

#endif
