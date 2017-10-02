#ifndef CepGen_StructureFunctions_Schaefer_h
#define CepGen_StructureFunctions_Schaefer_h

#include "StructureFunctions.h"

extern "C"
{
  extern void F2_fit_luxlike_( double& xbj, double& q2, double& F2, double& FL );
}

namespace CepGen
{
  namespace SF
  {
    StructureFunctions Schaefer( double q2, double xbj );
  }
}

#endif
