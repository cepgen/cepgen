#ifndef CepGen_Processes_ProcessesHandler_h
#define CepGen_Processes_ProcessesHandler_h

#include "CepGen/Core/ModuleFactory.h"
#include "CepGen/Processes/GenericProcess.h"

#define BUILDERNM( obj ) obj ## Builder
#define STRINGIFY( name ) #name
#define REGISTER_PROCESS( name, obj ) \
  struct BUILDERNM( name ) { \
    BUILDERNM( name )() { cepgen::ProcessesHandler::get().registerModule( STRINGIFY( name ), new obj ); } }; \
  static BUILDERNM( name ) g ## name;

namespace cepgen
{
  typedef ModuleFactory<cepgen::proc::GenericProcess> ProcessesHandler;
}

#endif
