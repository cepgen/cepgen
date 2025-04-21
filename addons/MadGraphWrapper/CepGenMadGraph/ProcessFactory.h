/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2025  Laurent Forthomme
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

#ifndef CepGenMadGraph_ProcessFactory_h
#define CepGenMadGraph_ProcessFactory_h

#include "CepGen/Modules/ModuleFactory.h"

/// Add a MadGraph process definition to the factory
#define REGISTER_MG5AMC_PROCESS(name, obj)                                     \
  namespace cepgen::mg5amc {                                                   \
    struct BUILDER_NAME(obj) {                                                 \
      BUILDER_NAME(obj)() { ProcessFactory::get().registerModule<obj>(name); } \
    };                                                                         \
    static const BUILDER_NAME(obj) gMGProc##obj;                               \
  }                                                                            \
  static_assert(true, "")

namespace cepgen::mg5amc {
  class Process;
  /// A MadGraph process factory
  DEFINE_FACTORY(ProcessFactory, Process, "MadGraph process definition factory");
}  // namespace cepgen::mg5amc

#endif
