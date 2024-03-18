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

#ifndef CepGen_Modules_ProcessFactory_h
#define CepGen_Modules_ProcessFactory_h

#include "CepGen/Modules/ModuleFactory.h"

/// Add a generic process definition to the list of handled processes
#define REGISTER_PROCESS(name, obj)                                           \
  namespace cepgen {                                                          \
    namespace proc {                                                          \
      struct BUILDERNM(obj) {                                                 \
        BUILDERNM(obj)() { ProcessFactory::get().registerModule<obj>(name); } \
      };                                                                      \
      static const BUILDERNM(obj) gProc##obj;                                 \
    }                                                                         \
  }                                                                           \
  static_assert(true, "")
/// Declare a Fortran process function name
#define DECLARE_FORTRAN_FUNCTION(f77_func) \
  extern "C" {                             \
  extern double f77_func##_();             \
  }                                        \
  static_assert(true, "")
/// Name of the process (Fortran scope)
#define PROCESS_F77_NAME(name) F77_##name
/// Convert a token into a string
#define STRINGIFY(name) #name
/// Add the Fortran process definition to the list of handled processes
#define REGISTER_FORTRAN_PROCESS(name, descr, f77_func)                           \
  struct PROCESS_F77_NAME(name) : public cepgen::proc::FortranFactorisedProcess { \
    PROCESS_F77_NAME(name)                                                        \
    (const cepgen::ParametersList& params = cepgen::ParametersList())             \
        : cepgen::proc::FortranFactorisedProcess(params, f77_func##_) {           \
      cepgen::proc::FortranFactorisedProcess::kProcParameters = params;           \
    }                                                                             \
    static cepgen::ParametersDescription description() {                          \
      auto desc = cepgen::proc::FortranFactorisedProcess::description();          \
      desc.setDescription(descr);                                                 \
      return desc;                                                                \
    }                                                                             \
  };                                                                              \
  REGISTER_PROCESS(STRINGIFY(name), F77_##name)

namespace cepgen {
  namespace proc {
    class Process;
  }
  /// A processes factory
  DEFINE_FACTORY(std::string, ProcessFactory, proc::Process, "Physics processes factory");
}  // namespace cepgen

#endif
