#ifndef CepGen_Processes_FortranProcesses_h
#define CepGen_Processes_FortranProcesses_h

#include "CepGen/Processes/FortranKTProcess.h"
#include "CepGen/Processes/ProcessesHandler.h"

#include "CepGen/Core/ParametersList.h"

#define DECLARE_FORTRAN_SUBROUTINE( method ) \
  extern "C" { extern void method ## _( double& ); }
#define PROCESS_F77_NAME( name ) F77_ ## name
#define REGISTER_FORTRAN_PROCESS( name, method, description ) \
  struct PROCESS_F77_NAME( name ) : public cepgen::proc::FortranKTProcess { \
    PROCESS_F77_NAME( name )() : cepgen::proc::FortranKTProcess( cepgen::ParametersList(), STRINGIFY( name ), description, method ## _ ) {} }; \
  REGISTER_PROCESS( name, PROCESS_F77_NAME( name ) )

#endif
