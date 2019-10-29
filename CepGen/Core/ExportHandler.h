#ifndef CepGen_IO_ExportHandler_h
#define CepGen_IO_ExportHandler_h

#include "CepGen/Core/ModuleFactory.h"
#include "CepGen/Core/GenericExportHandler.h"

#define REGISTER_IO_MODULE( name, obj ) \
  namespace cepgen { namespace io { \
    struct BUILDERNM( obj ) { \
      BUILDERNM( obj )() { ExportHandler::get().registerModule<obj>( name ); } }; \
    static BUILDERNM( obj ) g ## obj; \
  } }

namespace cepgen
{
  namespace io
  {
    /// An output modules factory
    typedef ModuleFactory<GenericExportHandler> ExportHandler;
  }
}

#endif

