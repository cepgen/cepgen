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

#ifndef CepGen_Modules_ProcessFactory_h
#define CepGen_Modules_ProcessFactory_h

#include "CepGen/Modules/ModuleFactory.h"

/** \file */

/// Add a generic process definition to the list of handled processes
#define REGISTER_PROCESS(name, obj)                                                                     \
  namespace cepgen {                                                                                    \
    namespace proc {                                                                                    \
      struct BUILDERNM(obj) {                                                                           \
        BUILDERNM(obj)() { ProcessFactory::get().registerModule<obj>(name, cepgen::ParametersList()); } \
      };                                                                                                \
      static const BUILDERNM(obj) gProc##obj;                                                           \
    }                                                                                                   \
  }
/// Declare a Fortran process function name
#define DECLARE_FORTRAN_FUNCTION(f77_func) \
  extern "C" {                             \
  extern double f77_func##_();             \
  }
#define PROCESS_F77_NAME(name) F77_##name
#define STRINGIFY(name) #name
/// Add the Fortran process definition to the list of handled processes
#define REGISTER_FORTRAN_PROCESS(name, descr, f77_func)                   \
  struct PROCESS_F77_NAME(name) : public cepgen::proc::FortranKTProcess { \
    PROCESS_F77_NAME(name)                                                \
    (const cepgen::ParametersList& params = cepgen::ParametersList())     \
        : cepgen::proc::FortranKTProcess(params, f77_func##_) {           \
      cepgen::proc::FortranKTProcess::kProcParameters = params;           \
    }                                                                     \
    static std::string description() { return descr; }                    \
  };                                                                      \
  REGISTER_PROCESS(STRINGIFY(name), F77_##name)

namespace cepgen {
  namespace proc {
    class Process;
    /// A processes factory
    DEFINE_FACTORY_STR(ProcessFactory, Process, "Physics processes factory");
  }  // namespace proc
}  // namespace cepgen

#endif
