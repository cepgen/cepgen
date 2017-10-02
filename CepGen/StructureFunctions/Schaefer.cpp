#include "Schaefer.h"

namespace CepGen
{
  namespace SF
  {
    StructureFunctions
    Schaefer( double q2, double xbj )
    {
      StructureFunctions luxlike;
#ifndef SchaeferF2
      FatalError( "LUXlike structure functions cannot be computed "
                  "as the Fortran subroutine is not linked to this instance!" );
#else
      F2_fit_luxlike_( xbj, q2, luxlike.F2, luxlike.FL );
#endif
      return luxlike;
    }
  }
}
