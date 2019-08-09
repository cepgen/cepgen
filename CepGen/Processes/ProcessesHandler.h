#ifndef CepGen_Processes_ProcessesHandler_h
#define CepGen_Processes_ProcessesHandler_h

#include "CepGen/Core/GenericProcess.h"
#include "CepGen/Core/ModuleFactory.h"
#include "CepGen/Core/ParametersList.h"

/** \file */

/// Add a generic process definition to the list of handled processes
#define REGISTER_PROCESS( name, obj ) \
  namespace cepgen { namespace proc { \
    struct BUILDERNM( obj ) { \
      BUILDERNM( obj )() { ProcessesHandler::get().registerModule<obj>( name ); } }; \
    static BUILDERNM( obj ) g ## obj; \
  } }
/// Declare a Fortran process function name
#define DECLARE_FORTRAN_FUNCTION( f77_func ) \
  extern "C" { extern double f77_func ## _(); }
#define PROCESS_F77_NAME( name ) F77_ ## name
#define STRINGIFY( name ) #name
/// Add the Fortran process definition to the list of handled processes
#define REGISTER_FORTRAN_PROCESS( name, f77_func, description ) \
  struct PROCESS_F77_NAME( name ) : public cepgen::proc::FortranKTProcess { \
    PROCESS_F77_NAME( name )( const cepgen::ParametersList& params = cepgen::ParametersList() ) : \
      cepgen::proc::FortranKTProcess( params, STRINGIFY( name ), description, f77_func ## _ ) { \
      cepgen::proc::FortranKTProcess::kProcParameters = params; \
    } }; \
  REGISTER_PROCESS( STRINGIFY( name ), F77_ ## name )

namespace cepgen
{
  namespace proc
  {
    /// A processes factory
    typedef ModuleFactory<GenericProcess> ProcessesHandler;
  }
}

#endif

