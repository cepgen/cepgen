/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2025  Laurent Forthomme
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

#ifndef CepGen_Modules_BaseIntegratorFactory_h
#define CepGen_Modules_BaseIntegratorFactory_h

#include "CepGen/Modules/ModuleFactory.h"

/// Add a generic integrator object builder definition
#define REGISTER_BASE_INTEGRATOR(name, obj)                                                \
  namespace cepgen {                                                                       \
    struct BUILDERNM(obj) {                                                                \
      BUILDERNM(obj)() { cepgen::BaseIntegratorFactory::get().registerModule<obj>(name); } \
    };                                                                                     \
    static const BUILDERNM(obj) gAnaInteg##obj;                                            \
  }                                                                                        \
  static_assert(true, "")

namespace cepgen {
  class BaseIntegrator;
  /// An integrator objects factory
  DEFINE_FACTORY(BaseIntegratorFactory, BaseIntegrator, "Integrators factory");
}  // namespace cepgen

#endif
