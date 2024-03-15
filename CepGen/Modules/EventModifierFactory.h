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
/** \file */

#ifndef CepGen_Modules_EventModifierFactory_h
#define CepGen_Modules_EventModifierFactory_h

#include "CepGen/Modules/ModuleFactory.h"

/// Add a generic event modification module definition to the factory
#define REGISTER_MODIFIER(name, obj)                                              \
  namespace cepgen {                                                              \
    struct BUILDERNM(obj) {                                                       \
      BUILDERNM(obj)() { EventModifierFactory::get().registerModule<obj>(name); } \
    };                                                                            \
    static const BUILDERNM(obj) gEveMod##obj;                                     \
  }                                                                               \
  static_assert(true, "")

namespace cepgen {
  class EventModifier;
  /// A event modifier algorithms factory
  DEFINE_FACTORY(std::string, EventModifierFactory, EventModifier, "Event modifiers factory");
}  // namespace cepgen

#endif
