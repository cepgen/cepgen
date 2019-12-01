#ifndef CepGen_Core_ExportHandler_h
#define CepGen_Core_ExportHandler_h

#include "CepGen/Core/ModuleFactory.h"
#include "CepGen/Modules/ExportModule.h"

#define REGISTER_IO_MODULE( name, obj ) \
  namespace cepgen { namespace io { \
    struct BUILDERNM( obj ) { \
      BUILDERNM( obj )() { ExportModuleHandler::get().registerModule<obj>( name ); } }; \
    static BUILDERNM( obj ) g ## obj; \
  } }

namespace cepgen
{
  namespace io
  {
    /// An output modules factory
    typedef ModuleFactory<ExportModule> ExportModuleHandler;
  }
}

#endif

