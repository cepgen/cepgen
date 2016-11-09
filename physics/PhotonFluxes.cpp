#include "PhotonFluxes.h"

double
PhotonFluxes::ProtonElastic( double x_, double kt2_ )
{
  double f_ela;

  const double mp = Particle::GetMassFromPDGId(Particle::Proton),
               mp2 = mp*mp;

  const double Q2_ela = ( kt2_+x_*x_*mp2 )/( 1.-x_ );
  const FormFactors ela = ElasticFormFactors( Q2_ela, mp2 );
  
  const double ela1 = pow( kt2_/( kt2_+x_*x_*mp2 ), 2 );
  const double ela2 = ela.FE;
  //const double ela3 = 1.-(Q2_ela-kt2_)/Q2_ela;
  //const double ela3 = 1.-pow(x_, 2)*mp2/Q2_ela/(1.-x_);
  //f_ela = alpha_em/Constants::Pi*(1.-x_+pow(x_, 2)/4.)*ela1*ela2*ela3/kt2_;
  f_ela = Constants::AlphaEM/Constants::Pi*ela1*ela2/Q2_ela;
  //f_ela = Constants::AlphaEM/Constants::Pi*( ( 1.-x_ )*ela1*ela2*ela3 + x_*x_/2.*G_M*G_M )/kt2_;

  return f_ela;
}

#ifdef GRVPDF

double
PhotonFluxes::ProtonInelastic( double x_, double kt2_, double mx_ )
{
  double f_ine;

  const double mx2 = mx_*mx_,
               mp = Particle::GetMassFromPDGId( Particle::Proton ),
               mp2 = mp*mp;

  const double Q02 = 0.8; // introduced to shift the Q2 scale
  double term1, term2;
  double f_aux;

  // F2 structure function
  const double Q2min = 1. / ( 1.-x_ )*( x_*( mx2-mp2 ) + x_*x_*mp2 ),
               Q2 = kt2_ / ( 1.-x_ ) + Q2min;
  float x_Bjorken = Q2 / ( Q2+mx2-mp2 );

  float mu2 = Q2+Q02; // scale is shifted

  float xuv, xdv, xus, xds, xss, xg;
  grv95lo_( x_Bjorken, mu2, xuv, xdv, xus, xds, xss, xg );

  const double F2_aux = 4./9.*( xuv + 2.*xus )
                      + 1./9.*( xdv + 2.*xds )
                      + 1./9.*2.*xss;

  // F2 corrected for low Q^2 behaviour
  const double F2_corr = Q2 / ( Q2+Q02 ) * F2_aux;

  ///////term1 = pow(1.- x_/2.*(mx2-mp2+Q2)/Q2, 2);
  //term1 = (1.-x_*(mx2-mp2+Q2)/Q2);
  term1 = ( 1. - ( Q2-kt2_ ) / Q2 );
  //term1 = (1.-Q2min/Q2);
  //term1 = 1.;
  term2 = pow( kt2_ / ( kt2_+x_*(mx2-mp2)+x_*x_*mp2 ), 2 );

  f_aux = F2_corr/( mx2+Q2-mp2 )*term1*term2;

  f_ine = Constants::AlphaEM/Constants::Pi*( 1.-x_ )*f_aux/kt2_;

  return f_ine;
}

#endif
