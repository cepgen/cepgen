#ifndef CepGen_Processes_FortranProcesses_h
#define CepGen_Processes_FortranProcesses_h

#include "CepGen/Processes/FortranKTProcess.h"

namespace CepGen
{
  namespace Process
  {
    struct FortranProcess
    {
      const char* name;
      void ( *method )( double& );
      const char* description;
    };
    static std::vector<FortranProcess> kFortranProcesses;
  }
}
void generateCepGenFortranProcesses();

#define F77_NAME( name ) name ## _
#define BEGIN_FORTRAN_PROCESSES_ENUM void generateCepGenFortranProcesses() {
#define REGISTER_FORTRAN_PROCESS( name, method, description ) \
  extern void F77_NAME( method )( double& ); \
  CepGen::Process::kFortranProcesses.emplace_back( CepGen::Process::FortranProcess{ name, F77_NAME( method ), description } );
#define END_FORTRAN_PROCESSES_ENUM }

#endif

