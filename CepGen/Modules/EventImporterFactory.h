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

#ifndef CepGen_Modules_EventImporterFactory_h
#define CepGen_Modules_EventImporterFactory_h

#include "CepGen/Modules/ModuleFactory.h"

/// Add an event import module definition to the factory
#define REGISTER_EVENT_IMPORTER(name, obj)                                           \
  namespace cepgen {                                                                 \
    struct BUILDER_NAME(obj) {                                                       \
      BUILDER_NAME(obj)() { EventImporterFactory::get().registerModule<obj>(name); } \
    };                                                                               \
    static const BUILDER_NAME(obj) gIOImporter##obj;                                 \
  }                                                                                  \
  static_assert(true, "")

namespace cepgen {
  class EventImporter;
  /// An event import algorithm factory
  DEFINE_FACTORY(EventImporterFactory, EventImporter, "Event importers factory");
}  // namespace cepgen

#endif
