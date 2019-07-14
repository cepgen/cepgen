#ifndef CepGen_Processes_ProcessesHandler_h
#define CepGen_Processes_ProcessesHandler_h

#include "CepGen/Core/ModuleFactory.h"
#include "CepGen/Core/ParametersList.h"

#include "CepGen/Processes/GenericProcess.h"
#include "CepGen/Processes/FortranKTProcess.h"

/** \file */

/// Add a generic process definition to the list of handled processes
#define REGISTER_PROCESS( name, obj ) \
  namespace cepgen { namespace proc { \
    struct BUILDERNM( name ) { \
      BUILDERNM( name )() { ProcessesHandler::get().registerModule<obj>( STRINGIFY( name ) ); } }; \
    static BUILDERNM( name ) g ## name; \
  } }
/// Declare a Fortran process function name
#define DECLARE_FORTRAN_FUNCTION( fort_func ) \
  extern "C" { extern double fort_func ## _(); }
#define PROCESS_F77_NAME( name ) F77_ ## name
/// Add the Fortran process definition to the list of handled processes
#define REGISTER_FORTRAN_PROCESS( name, fort_func, description ) \
  struct PROCESS_F77_NAME( name ) : public cepgen::proc::FortranKTProcess { \
    PROCESS_F77_NAME( name )( const cepgen::ParametersList& params = cepgen::ParametersList() ) : \
      cepgen::proc::FortranKTProcess( params, STRINGIFY( name ), description, fort_func ## _ ) { \
      cepgen::proc::FortranKTProcess::kProcParameters[STRINGIFY( fort_func )] = params; \
    } }; \
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

