#include "SzczurekUleshchenko.h"

namespace CepGen
{
  namespace SF
  {
    StructureFunctions
    SzczurekUleshchenko( double q2, double xbj )
    {
#ifndef GRVPDF
      FatalError( "Szczurek-Uleshchenko structure functions cannot be computed"
                  " as GRV PDF set is not linked to this instance!" );
#else
      const float q02 = 0.8;
      float amu2 = q2+q02; // shift the overall scale
      float xuv, xdv, xus, xds, xss, xg;
      float xbj_arg = xbj;

      grv95lo_( xbj_arg, amu2, xuv, xdv, xus, xds, xss, xg );

      DebuggingInsideLoop( Form( "Form factor content at xB = %e (scale = %f GeV^2):\n\t"
                                 "  valence quarks: u / d     = %e / %e\n\t"
                                 "  sea quarks:     u / d / s = %e / %e / %e\n\t"
                                 "  gluons:                   = %e",
                                 xbj, amu2, xuv, xdv, xus, xds, xss, xg ) );

      // standard partonic structure function
      const double F2_aux = 4./9.*( xuv + 2.*xus )
                          + 1./9.*( xdv + 2.*xds )
                          + 1./9.*(       2.*xss );

      // F2 corrected for low Q^2 behaviour
      const double F2_corr = F2_aux * q2 / amu2,
                   F1 = 0.5*F2_corr/xbj; // Callan-Gross relation

      return StructureFunctions( F1, F2_corr );
#endif
    }
  }
}
