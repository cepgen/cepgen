#ifndef CepGen_Module_ExportModuleFactory_h
#define CepGen_Module_ExportModuleFactory_h

#include "CepGen/Modules/ModuleFactory.h"

/** \file */

/// Add a generic export module definition to the factory
#define REGISTER_IO_MODULE(name, obj)                                              \
  namespace cepgen {                                                               \
    namespace io {                                                                 \
      struct BUILDERNM(obj) {                                                      \
        BUILDERNM(obj)() { ExportModuleFactory::get().registerModule<obj>(name); } \
      };                                                                           \
      static BUILDERNM(obj) gIO##obj;                                              \
    }                                                                              \
  }

namespace cepgen {
  namespace io {
    class ExportModule;
    /// An output modules factory
    typedef ModuleFactory<ExportModule> ExportModuleFactory;
  }  // namespace io
}  // namespace cepgen

#endif
