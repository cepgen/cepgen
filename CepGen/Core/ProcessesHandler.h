#ifndef CepGen_Core_ProcessesHandler_h
#define CepGen_Core_ProcessesHandler_h

#include "CepGen/Core/ModuleFactory.h"
#include "CepGen/Processes/GenericProcess.h"
#include "CepGen/Processes/FortranKTProcess.h"
#include "CepGen/Core/ParametersList.h"

#define BUILDERNM( obj ) obj ## Builder
#define STRINGIFY( name ) #name
#define REGISTER_PROCESS( name, obj ) \
  struct BUILDERNM( name ) { \
    BUILDERNM( name )() { cepgen::proc::ProcessesHandler::get().registerModule<obj>( STRINGIFY( name ) ); } }; \
  static BUILDERNM( name ) g ## name;
#define DECLARE_FORTRAN_SUBROUTINE( method ) \
  extern "C" { extern void method ## _( double& ); }
#define PROCESS_F77_NAME( name ) F77_ ## name
#define REGISTER_FORTRAN_PROCESS( name, method, description ) \
  struct PROCESS_F77_NAME( name ) : public cepgen::proc::FortranKTProcess { \
    PROCESS_F77_NAME( name )( const cepgen::ParametersList& params = cepgen::ParametersList() ) : cepgen::proc::FortranKTProcess( params, STRINGIFY( name ), description, method ## _ ) {} }; \
  REGISTER_PROCESS( name, PROCESS_F77_NAME( name ) )

namespace cepgen
{
  namespace proc
  {
    typedef ModuleFactory<GenericProcess> ProcessesHandler;
  }
}

#endif

