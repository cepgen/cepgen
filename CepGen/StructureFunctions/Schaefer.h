#ifndef CepGen_StructureFunctions_Schaefer_h
#define CepGen_StructureFunctions_Schaefer_h

#include "StructureFunctions.h"

extern "C"
{
  extern void f2_fit_luxlike_( double& xbj, double& q2, double& F2, double& FL );
}

namespace CepGen
{
  namespace SF
  {
    class Schaefer
    {
      public:
        Schaefer() {}
        StructureFunctions operator()( double q2, double xbj ) const;
    };
  }
}

#endif
