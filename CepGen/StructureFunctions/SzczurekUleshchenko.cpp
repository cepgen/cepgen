#include "CepGen/StructureFunctions/SzczurekUleshchenko.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"

namespace cepgen
{
  namespace strfun
  {
    SzczurekUleshchenko::SzczurekUleshchenko() :
      Parameterisation( Type::SzczurekUleshchenko ), F1( 0. )
    {}

    SzczurekUleshchenko&
    SzczurekUleshchenko::operator()( double xbj, double q2 )
    {
#ifndef GRVPDF
      throw CG_FATAL( "SzczurekUleshchenko" )
        << "Szczurek-Uleshchenko structure functions cannot be computed"
        << " as GRV PDF set is not linked to this instance!";
#else
      std::pair<double,double> nv = { xbj, q2 };
      if ( nv == old_vals_ )
        return *this;
      old_vals_ = nv;

      const float q02 = 0.8;
      float amu2 = q2+q02; // shift the overall scale
      float xuv, xdv, xus, xds, xss, xg;
      float xbj_arg = xbj;

      grv95lo_( xbj_arg, amu2, xuv, xdv, xus, xds, xss, xg );

      CG_DEBUG_LOOP( "SzczurekUleshchenko" )
        << "Form factor content at xB = " << xbj << " (scale = " << amu2 << " GeV^2):\n\t"
        << "  valence quarks: u / d     = " << xuv << " / " << xdv << "\n\t"
        << "  sea quarks:     u / d / s = " << xus << " / " << xds << " / " << xss << "\n\t"
        << "  gluons:                   = " << xg;

      // standard partonic structure function
      const double F2_aux = 4./9.*( xuv + 2.*xus )
                          + 1./9.*( xdv + 2.*xds )
                          + 1./9.*(       2.*xss );

      F2 = F2_aux * q2 / amu2; // F2 corrected for low Q^2 behaviour

      return *this;
#endif
    }
  }
}
