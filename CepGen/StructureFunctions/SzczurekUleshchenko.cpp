#include "SzczurekUleshchenko.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"

namespace CepGen
{
  namespace SF
  {
    SzczurekUleshchenko
    SzczurekUleshchenko::operator()( double q2, double xbj, const SigmaRatio& ratio ) const
    {
#ifndef GRVPDF
      throw FatalError( "SzczurekUleshchenko" )
        << "Szczurek-Uleshchenko structure functions cannot be computed"
        << " as GRV PDF set is not linked to this instance!";
#else
      const float q02 = 0.8;
      float amu2 = q2+q02; // shift the overall scale
      float xuv, xdv, xus, xds, xss, xg;
      float xbj_arg = xbj;

      grv95lo_( xbj_arg, amu2, xuv, xdv, xus, xds, xss, xg );

      DebuggingInsideLoop( "SzczurekUleshchenko" )
        << "Form factor content at xB = " << xbj << " (scale = " << amu2 << " GeV^2):\n\t"
        << "  valence quarks: u / d     = " << xuv << " / " << xdv << "\n\t"
        << "  sea quarks:     u / d / s = " << xus << " / " << xds << " / " << xss << "\n\t"
        << "  gluons:                   = " << xg;

      // standard partonic structure function
      const double F2_aux = 4./9.*( xuv + 2.*xus )
                          + 1./9.*( xdv + 2.*xds )
                          + 1./9.*(       2.*xss );

      SzczurekUleshchenko su;
      su.F2 = F2_aux * q2 / amu2; // F2 corrected for low Q^2 behaviour
      double r_err = 0.;
      su.computeFL( q2, xbj, ratio( q2, xbj, r_err ) );

      return su;
#endif
    }
  }
}
