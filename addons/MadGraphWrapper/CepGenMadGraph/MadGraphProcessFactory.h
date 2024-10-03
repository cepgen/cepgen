/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2024  Laurent Forthomme
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

#ifndef CepGenMadGraph_MadGraphProcessFactory_h
#define CepGenMadGraph_MadGraphProcessFactory_h

#include "CepGen/Modules/ModuleFactory.h"

/// Add a MadGraph process definition to the factory
#define REGISTER_MG5AMC_PROCESS(name, obj)                                          \
  namespace cepgen {                                                                \
    struct BUILDERNM(obj) {                                                         \
      BUILDERNM(obj)() { MadGraphProcessFactory::get().registerModule<obj>(name); } \
    };                                                                              \
    static const BUILDERNM(obj) gMGProc##obj;                                       \
  }                                                                                 \
  static_assert(true, "")

namespace cepgen {
  class MadGraphProcess;
  /// A MadGraph process factory
  DEFINE_FACTORY(MadGraphProcessFactory, MadGraphProcess, "MadGraph process definition factory");
}  // namespace cepgen

#endif
