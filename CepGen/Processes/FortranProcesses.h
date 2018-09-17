#ifndef CepGen_Processes_FortranProcesses_h
#define CepGen_Processes_FortranProcesses_h

#include "CepGen/Processes/FortranKTProcess.h"

#ifdef __cplusplus__
extern "C" {
#endif
  extern void nucl_to_ff_( double& );
#ifdef __cplusplus__
}
#endif

namespace CepGen
{
  namespace Process
  {
    static std::map<const char*,void(*)( double& )> kFortranProcesses;
  }
}

#endif

