#include "CepGen/StructureFunctions/StructureFunctions.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"

extern "C"
{
  extern void grv95lo_( float&, float&, float&, float&, float&, float&, float&, float& );
}

namespace cepgen
{
  namespace strfun
  {
    /// Szcurek and Uleshchenko modelling of \f$F_{1,2}\f$ \cite Szczurek:1999wp
    class SzczurekUleshchenko : public Parameterisation
    {
      public:
        SzczurekUleshchenko( const ParametersList& params = ParametersList() );
        SzczurekUleshchenko& operator()( double xbj, double q2 ) override;

      private:
        /// \f$Q^2\f$ scale shift
        const float q2_shift_;

      public:
        double F1;
    };

    SzczurekUleshchenko::SzczurekUleshchenko( const ParametersList& params ) :
      Parameterisation( params ),
      q2_shift_( params.get<double>( "q2shift", 0.8 ) ),
      F1( 0. )
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

      float amu2 = q2+q2_shift_; // shift the overall scale
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

REGISTER_STRFUN( SzczurekUleshchenko, strfun::SzczurekUleshchenko )
