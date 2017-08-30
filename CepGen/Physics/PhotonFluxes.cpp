#include "PhotonFluxes.h"

namespace CepGen
{
  namespace Fluxes
  {
    namespace Photon
    {
      double
      ProtonElastic( double x, double kt2 )
      {
        const double mp = Particle::massFromPDGId( Particle::Proton ), mp2 = mp*mp;

        const double Q2_ela = ( kt2+x*x*mp2 )/( 1.-x );
        const FormFactors ela = ElasticFormFactors( Q2_ela, mp2 );

        const double ela1 = pow( kt2/( kt2+x*x*mp2 ), 2 );
        const double ela2 = ela.FE;
        //const double ela3 = 1.-(Q2_ela-kt2)/Q2_ela;
        //const double ela3 = 1.-x*x*mp2/Q2_ela/(1.-x);
        //return Constants::alpha_em/M_PI*(1.-x+x*x/4.)*ela1*ela2*ela3/kt2;
        return Constants::alphaEM/M_PI*ela1*ela2/Q2_ela;
        //return Constants::alphaEM/M_PI*( ( 1.-x )*ela1*ela2*ela3 + x*x/2.*G_M*G_M )/kt2;
      }


      double
      ProtonInelastic( double x, double kt2, double mx )
      {
#ifndef GRVPDF
        FatalError( "Inelastic flux cannot be computed as GRV PDF set is not linked to this instance!" );
#else
        const double mx2 = mx*mx,
                     mp = Particle::massFromPDGId( Particle::Proton ),
                     mp2 = mp*mp;

        const double Q02 = 0.8; // introduced to shift the Q2 scale
        double term1, term2;
        double f_aux;

        // F2 structure function
        const double Q2min = 1. / ( 1.-x )*( x*( mx2-mp2 ) + x*x*mp2 ),
                     Q2 = kt2 / ( 1.-x ) + Q2min;
        float x_Bjorken = Q2 / ( Q2+mx2-mp2 );

        float mu2 = Q2+Q02; // scale is shifted

        float xuv, xdv, xus, xds, xss, xg;
        grv95lo_( x_Bjorken, mu2, xuv, xdv, xus, xds, xss, xg );

        const double F2_aux = 4./9.*( xuv + 2.*xus )
                            + 1./9.*( xdv + 2.*xds )
                            + 1./9.*2.*xss;

        // F2 corrected for low Q^2 behaviour
        const double F2_corr = Q2 / ( Q2+Q02 ) * F2_aux;

        ///////term1 = pow(1.- x/2.*(mx2-mp2+Q2)/Q2, 2);
        //term1 = (1.-x*(mx2-mp2+Q2)/Q2);
        term1 = ( 1. - ( Q2-kt2 ) / Q2 );
        //term1 = (1.-Q2min/Q2);
        //term1 = 1.;
        term2 = pow( kt2 / ( kt2+x*(mx2-mp2)+x*x*mp2 ), 2 );

        f_aux = F2_corr/( mx2+Q2-mp2 )*term1*term2;

        return Constants::alphaEM/M_PI*( 1.-x )*f_aux/kt2;
#endif
      }
    }
  }
}

