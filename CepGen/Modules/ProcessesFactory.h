#ifndef CepGen_Processes_ProcessesFactory_h
#define CepGen_Processes_ProcessesFactory_h

#include "CepGen/Core/ModuleFactory.h"

/** \file */

/// Add a generic process definition to the list of handled processes
#define REGISTER_PROCESS( name, obj ) \
  namespace cepgen { namespace proc { \
    struct BUILDERNM( obj ) { \
      BUILDERNM( obj )() { ProcessesFactory::get().registerModule<obj>( name, \
        cepgen::ParametersList() ); } }; \
    static BUILDERNM( obj ) gProc ## obj; \
  } }
/// Declare a Fortran process function name
#define DECLARE_FORTRAN_FUNCTION( f77_func ) \
  extern "C" { extern double f77_func ## _(); }
#define PROCESS_F77_NAME( name ) F77_ ## name
#define STRINGIFY( name ) #name
/// Add the Fortran process definition to the list of handled processes
#define REGISTER_FORTRAN_PROCESS( name, descr, f77_func ) \
  struct PROCESS_F77_NAME( name ) : public cepgen::proc::FortranKTProcess { \
    PROCESS_F77_NAME( name )( const cepgen::ParametersList& params = cepgen::ParametersList() ) : \
      cepgen::proc::FortranKTProcess( params, f77_func ## _ ) { \
      cepgen::proc::FortranKTProcess::kProcParameters = params; \
    } \
    static std::string description() { return descr; } \
  }; \
  REGISTER_PROCESS( STRINGIFY( name ), F77_ ## name )

namespace cepgen
{
  namespace proc
  {
    class Process;
    /// A processes factory
    typedef ModuleFactory<Process> ProcessesFactory;
  }
}

#endif

