/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2023  Laurent Forthomme
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

#ifndef CepGen_Modules_DerivatorFactory_h
#define CepGen_Modules_DerivatorFactory_h

#include "CepGen/Modules/ModuleFactory.h"

/** \file */

/// Add a generic derivator object builder definition
#define REGISTER_DERIVATOR(name, obj)                                           \
  namespace cepgen {                                                            \
    namespace utils {                                                           \
      struct BUILDERNM(obj) {                                                   \
        BUILDERNM(obj)() { DerivatorFactory::get().registerModule<obj>(name); } \
      };                                                                        \
      static const BUILDERNM(obj) gDeriv##obj;                                  \
    }                                                                           \
  }                                                                             \
  static_assert(true, "")

namespace cepgen {
  namespace utils {
    class Derivator;
  }
  /// An analytical integrator objects factory
  DEFINE_FACTORY_STR(DerivatorFactory, utils::Derivator, "Derivators factory");
}  // namespace cepgen

#endif
