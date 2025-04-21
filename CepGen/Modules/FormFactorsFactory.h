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

#ifndef CepGen_Modules_FormFactorsFactory_h
#define CepGen_Modules_FormFactorsFactory_h

#include "CepGen/Modules/ModuleFactory.h"

/// Add an electromagnetic form factors definition to the list of handled parameterisations
#define REGISTER_FORMFACTORS(name, obj)                                            \
  namespace cepgen::formfac {                                                      \
    struct BUILDER_NAME(obj) {                                                     \
      BUILDER_NAME(obj)() { FormFactorsFactory::get().registerModule<obj>(name); } \
    };                                                                             \
    static const BUILDER_NAME(obj) gFF##obj;                                       \
  }                                                                                \
  static_assert(true, "")

namespace cepgen::formfac {
  class Parameterisation;
  /// Standard dipole handler name
  static constexpr auto gFFStandardDipoleHandler = "StandardDipole";
}  // namespace cepgen::formfac

namespace cepgen {
  /// A form factors parameterisations factory
  DEFINE_FACTORY(FormFactorsFactory, formfac::Parameterisation, "Nucleon form factors factory");
}  // namespace cepgen

#endif
