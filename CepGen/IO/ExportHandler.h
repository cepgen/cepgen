#ifndef CepGen_IO_ExportHandler_h
#define CepGen_IO_ExportHandler_h

#include "CepGen/Core/ModuleFactory.h"
#include "CepGen/IO/GenericExportHandler.h"

#define REGISTER_IO_MODULE( name, obj ) \
  namespace cepgen { namespace io { \
    struct BUILDERNM( name ) { \
      BUILDERNM( name )() { ExportHandler::get().registerModule<obj>( STRINGIFY( name ) ); } }; \
    static BUILDERNM( name ) g ## name; \
  } }

namespace cepgen
{
  namespace io
  {
    /// A hadroniser modules factory
    typedef ModuleFactory<GenericExportHandler> ExportHandler;
  }
}

#endif

