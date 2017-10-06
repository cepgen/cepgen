#include "Schaefer.h"

namespace CepGen
{
  namespace SF
  {
    StructureFunctions
    Schaefer::operator()( double q2, double xbj ) const
    {
      StructureFunctions luxlike;
#ifndef SchaeferF2
      FatalError( "LUXlike structure functions cannot be computed "
                  "as the Fortran subroutine is not linked to this instance!" );
#else
      f2_fit_luxlike_( xbj, q2, luxlike.F2, luxlike.FL );
#endif
      return luxlike;
    }
  }
}
