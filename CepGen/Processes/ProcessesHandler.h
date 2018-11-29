#ifndef CepGen_Processes_ProcessesHandler_h
#define CepGen_Processes_ProcessesHandler_h

#include "CepGen/Core/ModuleFactory.h"
#include "CepGen/Core/ParametersList.h"

#include "CepGen/Processes/GenericProcess.h"
#include "CepGen/Processes/FortranKTProcess.h"

#define BUILDERNM( obj ) obj ## Builder
#define STRINGIFY( name ) #name
#define REGISTER_PROCESS( name, obj ) \
  namespace cepgen { namespace proc { \
    struct BUILDERNM( name ) { \
      BUILDERNM( name )() { ProcessesHandler::get().registerModule<obj>( STRINGIFY( name ) ); } }; \
    static BUILDERNM( name ) g ## name; \
  } }
#define DECLARE_FORTRAN_FUNCTION( method ) \
  extern "C" { extern double method ## _(); }
#define PROCESS_F77_NAME( name ) F77_ ## name
#define REGISTER_FORTRAN_PROCESS( name, method, description ) \
  struct PROCESS_F77_NAME( name ) : public cepgen::proc::FortranKTProcess { \
    PROCESS_F77_NAME( name )( const cepgen::ParametersList& params = cepgen::ParametersList() ) : cepgen::proc::FortranKTProcess( params, STRINGIFY( name ), description, method ## _ ) {} }; \
  REGISTER_PROCESS( name, PROCESS_F77_NAME( name ) )

namespace cepgen
{
  namespace proc
  {
    /// A processes factory
    typedef ModuleFactory<GenericProcess> ProcessesHandler;
  }
}

#endif

