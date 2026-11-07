/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2026  Laurent Forthomme
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

#ifndef CepGen_Modules_SurvivalFactorFactory_h
#define CepGen_Modules_SurvivalFactorFactory_h

#include "CepGen/Modules/ModuleFactory.h"

/// Add a gap survival factor modelling to the list of handled algorithms
#define REGISTER_SURVIVAL_FACTOR(name, obj)                                           \
  namespace cepgen {                                                                  \
    struct BUILDER_NAME(obj) {                                                        \
      BUILDER_NAME(obj)() { SurvivalFactorFactory::get().registerModule<obj>(name); } \
    };                                                                                \
    static const BUILDER_NAME(obj) gSurvFact##obj;                                    \
  }                                                                                   \
  static_assert(true, "")

namespace cepgen {
  class SurvivalFactor;
  /// A gap survival factor modellings factory
  DEFINE_FACTORY(SurvivalFactorFactory, SurvivalFactor, "Gap survival factor modelling factory");
}  // namespace cepgen

#endif
