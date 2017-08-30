#include "FormFactors.h"

namespace CepGen
{
  FormFactors
  TrivialFormFactors()
  {
    return FormFactors( 1.0, 1.0 );
  }

  FormFactors
  ElasticFormFactors( double q2, double mi2 )
  {
    const double GE = pow(1.+q2/0.71, -2.), GM = 2.79*GE;
    return FormFactors( ( 4.*mi2*GE*GE+q2*GM*GM ) / ( 4.*mi2 + q2 ), GM*GM );
  }

  FormFactors
  SuriYennieFormFactors( double q2, double mi2, double mf2 )
  {
    // values extracted from experimental fits
    const double cc1 = 0.86926, // 0.6303
                 cc2 = 2.23422, // 2.2049
                 dd1 = 0.12549, // 0.0468
                 cp = 0.96, // 1.23
                 bp = 0.63, // 0.61
                 rho = 0.585; // 1.05
    const double x = q2/(q2+mf2),
                 dm2 = mf2-mi2,
                 en = dm2+q2,
                 tau = -q2/4./mi2,
                 rhot = rho+q2;
    const double rho_norm = rho/rhot;

    FormFactors ff;
    ff.FM = ( -1./q2 ) * ( -cc1*rho_norm*rho_norm*dm2 - cc2*mi2*pow( 1.-x, 4 )/( x*( x*cp-2*bp )+1. ) );
    ff.FE = ( -tau*ff.FM + dd1*dm2*q2*rho_norm*pow( dm2/en, 2 )/( rhot*mi2 ) )/( 1. + en*en/( 4.*mi2*q2 ) );
    return ff;
  }

  FormFactors
  FioreBrasseFormFactors( double q2, double mi2, double mf2 )
  {
    const double k = 2.*sqrt( mi2 );
    // start by computing the proton structure function for this Q**2/mX couple
    double dummy, psfw1, psfw2;
    if ( !PSF( -q2, mf2, dummy, psfw1, psfw2 ) ) {
      InWarning( Form( "Fiore-Brasse form factors to be retrieved for an invalid MX value:\n\t"
                       "%.2e GeV, while allowed range is [1.07, 1.99] GeV", sqrt( mf2 ) ) );
      return FormFactors( 0.0, 0.0 );
    }

    return FormFactors( psfw2 / k, -psfw1*k / q2 );
  }

  FormFactors
  SzczurekUleshchenkoFormFactors( double q2, double mi2, double mf2 )
  {
    const double k = 2.*sqrt( mi2 ),
                 q2_0 = 0.8;

    float x = q2 / ( mf2+q2+mi2 ),
          amu2 = q2+q2_0; // shift the overall scale
    float xuv, xdv, xus, xds, xss, xg;

    grv95lo_( x, amu2, xuv, xdv, xus, xds, xss, xg );

    DebuggingInsideLoop( Form( "Form factor content at xB = %e (scale = %f GeV^2):\n\t"
                               "  valence quarks: u / d     = %e / %e\n\t"
                               "  sea quarks:     u / d / s = %e / %e / %e\n\t"
                               "  gluons:                   = %e",
                               x, amu2, xuv, xdv, xus, xds, xss, xg ) );

    const double F2_aux = 4./9.*( xuv+2.*xus )
                        + 1./9.*( xdv+2.*xds )
                        + 1./9.*2.*xss;
    /*const double F2_aux = 4./9.*( xuv + 2.*xus )
                        + 1./9.*( 0. + 2.*xds )
                        + 1./9.*2.*xss;*/

    // F2 corrected for low Q^2 behaviour
    const double F2_corr = q2/amu2*F2_aux,
                 F1 = F2_corr/( 2.*x ); // Callan-Gross relation

    /*const double F2_corr = Q2 / ( Q2+Q02 ) * F2_aux;
    ///////term1 = pow(1.- x_/2.*(mx2-mp2+Q2)/Q2, 2);
    //term1 = (1.-x_*(mx2-mp2+Q2)/Q2);
    term1 = ( 1. - ( Q2-kt2_ ) / Q2 );
    //term1 = (1.-Q2min/Q2);
    //term1 = 1.;
    term2 = pow( kt2_ / ( kt2_+x_*(mx2-mp2)+x_*x_*mp2 ), 2 );
    f_aux = F2_corr/( mx2+Q2-mp2 )*term1*term2;
    f_ine = Constants::alphaEM/M_PI*( 1.-x_ )*f_aux/kt2_;
    return f_ine;*/

    const double w2 = k*x/q2*F2_corr,
                 w1 = 2.*F1/k;

    return FormFactors( w2/k, -w1*k/q2 );
  }

  std::ostream&
  operator<<( std::ostream& os, const FormFactors& ff )
  {
    os << Form( "Form factors: electric: Fe = %.3e ; magnetic: Fm = %.3e", ff.FE, ff.FM ).c_str();
    return os;
  }
}

