#include "FormFactors.h"

bool
PSF( double q2, double mx2, double* sigma_t, double* w1, double* w2 )
{
  int nBin;
  double xBin, dx, nu2, logqq0, gd2;
  double sigLow, sigHigh;

  //const double m_min = Particle::GetMassFromPDGId(Particle::Proton)+0.135;
  const double m_proton = Particle::GetMassFromPDGId( Particle::Proton ),
               m2_proton = m_proton*m_proton,
               m_min = m_proton+Particle::GetMassFromPDGId( Particle::PiZero );

  const double mx = sqrt( mx2 );

  if ( mx<m_min or mx<1.99 ) {
    *sigma_t = *w1 = *w2 = 0.;
    return false;
  }

  // values of a, b, c provided from the fits on ep data and retrieved from
  // http://dx.doi.org/10.1016/0550-3213(76)90231-5 with 1.110 <= w2 <=1.990
  const double abrass[56] = { 5.045,5.126,5.390,5.621,5.913,5.955,6.139,6.178,6.125,5.999,
                              5.769,5.622,5.431,5.288,5.175,5.131,5.003,5.065,5.045,5.078,
                              5.145,5.156,5.234,5.298,5.371,5.457,5.543,5.519,5.465,5.384,
                              5.341,5.320,5.275,5.290,5.330,5.375,5.428,5.478,5.443,5.390,
                              5.333,5.296,5.223,5.159,5.146,5.143,5.125,5.158,5.159,5.178,
                              5.182,5.195,5.160,5.195,5.163,5.172 },
               bbrass[56] = { 0.798,1.052,1.213,1.334,1.397,1.727,1.750,1.878,1.887,1.927,
                              2.041,2.089,2.148,2.205,2.344,2.324,2.535,2.464,2.564,2.610,
                              2.609,2.678,2.771,2.890,2.982,3.157,3.183,3.315,3.375,3.450,
                              3.477,3.471,3.554,3.633,3.695,3.804,3.900,4.047,4.290,4.519,
                              4.709,4.757,4.840,5.017,5.015,5.129,5.285,5.322,5.545,5.623,
                              5.775,5.894,6.138,6.151,6.301,6.542 },
               cbrass[56] = { 0.043, 0.024, 0.000,-0.013,-0.023,-0.069,-0.060,-0.080,-0.065,-0.056,
                             -0.065,-0.056,-0.043,-0.034,-0.054,-0.018,-0.046,-0.015,-0.029,-0.048,
                             -0.032,-0.045,-0.084,-0.115,-0.105,-0.159,-0.164,-0.181,-0.203,-0.223,
                             -0.245,-0.254,-0.239,-0.302,-0.299,-0.318,-0.383,-0.393,-0.466,-0.588,
                             -0.622,-0.568,-0.574,-0.727,-0.665,-0.704,-0.856,-0.798,-1.048,-0.980,
                             -1.021,-1.092,-1.313,-1.341,-1.266,-1.473 };

  if ( mx<1.11 ) {
    nBin = 0;
    xBin = mx-m_min;
    dx = 1.11-m_min; // Delta w bin sizes
  }
  else if ( mx<1.77 ) { // w in [1.11, 1.77[
    dx = 0.015; // Delta w bin sizes
    nBin = ( mx-1.11 )/dx + 1;
    xBin = fmod( mx-1.11, dx );
  }
  else { // w in [1.77, 1.99[
    dx = 0.02; // Delta w bin sizes
    nBin = ( mx-1.77 )/dx + 45;
    xBin = fmod( mx-1.77, dx );
  }
  nu2 = pow( ( mx2-q2-m2_proton ) / ( 2.*m_proton ), 2 );
  logqq0 = log( ( nu2-q2 ) / pow( ( mx2-m2_proton ) / ( 2.*m_proton ), 2 ) ) / 2.;
  gd2 = pow( 1. / ( 1-q2 / .71 ), 4 ); // dipole form factor of the proton

  sigLow = (nBin==0) ? 0. :
    exp( abrass[nBin-1]+bbrass[nBin-1]*logqq0+cbrass[nBin-1]*pow( fabs( logqq0 ), 3 ) )*gd2;
  sigHigh =
    exp( abrass[nBin]  +bbrass[nBin]  *logqq0+cbrass[nBin]  *pow( fabs( logqq0 ), 3 ) )*gd2;

  *sigma_t = sigLow + xBin*( sigHigh-sigLow )/dx;
  *w1 = ( mx2-m2_proton )/( 8.*Constants::Pi*Constants::Pi*m_proton*Constants::AlphaEM )/Constants::GeV2toBarn*1.e6*(*sigma_t);
  *w2 = ( *w1 ) * q2/( q2-nu2 );

  return true;
}

FormFactors
TrivialFormFactors()
{
  FormFactors ff;
  ff.FE = 1.;
  ff.FM = 1.;
  return ff;
}

FormFactors
ElasticFormFactors( double q2, double mi2 )
{
  const double GE = pow(1.+q2/0.71, -2.), GM = 2.79*GE;
  FormFactors ff;
  ff.FE = ( 4.*mi2*GE*GE+q2*GM*GM ) / ( 4.*mi2 + q2 );
  ff.FM = GM*GM;
  return ff;
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
  double dummy, psfw1, psfw2; PSF( -q2, mf2, &dummy, &psfw1, &psfw2 );

  FormFactors ff;
  ff.FM =-psfw1*k / q2;
  ff.FE = psfw2 / k;
  return ff;
}

FormFactors
SzczurekUleschenkoFormFactors( double q2, double mi2, double mf2 )
{
  const double k = 2.*sqrt( mi2 ),
               q2_0 = 0.8;

  float x = q2 / ( mf2+q2+mi2 ),
        amu2 = q2+q2_0; // shift the overall scale
  float xuv, xdv, xus, xds, xss, xg;

  grv95lo_( x, amu2, xuv, xdv, xus, xds, xss, xg );

  const double F2_aux = 4./9.*(xuv+2.*xus)
                      + 1./9.*(xdv+2.*xds)
                      + 1./9.*2.*xss;

  // F2 corrected for low Q^2 behaviour
  const double F2_corr = q2/amu2*F2_aux,
               F1 = F2_corr/(2.*x); // Callan-Gross relation

  const double w2 = k*x/q2*F2_corr,
               w1 = 2.*F1/k;

  FormFactors ff;
  ff.FM =-w1*k / q2;
  ff.FE = w2 / k;
  return ff;
}

std::ostream&
operator<<( std::ostream& os, const FormFactors& ff )
{
  os << Form( "Form factors: electric: Fe = %.3e ; magnetic: Fm = %.3e", ff.FE, ff.FM ).c_str();
  return os;
}
