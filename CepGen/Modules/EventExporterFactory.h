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

#ifndef CepGen_Modules_EventExporterFactory_h
#define CepGen_Modules_EventExporterFactory_h

#include "CepGen/Modules/ModuleFactory.h"

/** \file */

/// Add a generic export module definition to the factory
#define REGISTER_EXPORTER(name, obj)                                              \
  namespace cepgen {                                                              \
    struct BUILDERNM(obj) {                                                       \
      BUILDERNM(obj)() { EventExporterFactory::get().registerModule<obj>(name); } \
    };                                                                            \
    static const BUILDERNM(obj) gIO##obj;                                         \
  }

namespace cepgen {
  class EventExporter;
  /// An output modules factory
  DEFINE_FACTORY_STR(EventExporterFactory, EventExporter, "Export modules factory");
}  // namespace cepgen

#endif
