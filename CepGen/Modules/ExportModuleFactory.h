#ifndef CepGen_Module_ExportModuleFactory_h
#define CepGen_Module_ExportModuleFactory_h

#include "CepGen/Core/ModuleFactory.h"

#define REGISTER_IO_MODULE( name, obj ) \
  namespace cepgen { namespace io { \
    struct BUILDERNM( obj ) { \
      BUILDERNM( obj )() { ExportModuleFactory::get().registerModule<obj>( name ); } }; \
    static BUILDERNM( obj ) gIO ## obj; \
  } }

namespace cepgen
{
  namespace io
  {
    class ExportModule;
    /// An output modules factory
    typedef ModuleFactory<ExportModule> ExportModuleFactory;
  }
}

#endif

