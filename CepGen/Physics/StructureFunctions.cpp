#include "StructureFunctions.h"

namespace CepGen
{
  std::ostream&
  operator<<( std::ostream& os, const StructureFunctions& sf )
  {
    switch ( sf ) {
      case Electron:            os << "electron"; break;
      case ElasticProton:       os << "elastic proton"; break;
      case SuriYennie:          os << "Suri-Yennie"; break;
      case SuriYennieLowQ2:     os << "Suri-Yennie;lowQ2"; break;
      case SzczurekUleshchenko: os << "Szczurek-Uleshchenko"; break;
      case FioreVal:            os << "Fiore;valence"; break;
      case FioreSea:            os << "Fiore;sea"; break;
      case Fiore:               os << "Fiore"; break;
    }
    return os;
  }

  /// Fiore-Brasse proton structure function
  StructureFunctions
  StructureFunctions::FioreBrasse( double q2, double mx2 )
  {
    sigma_t = w1 = w2 = 0.;

    //const double m_min = Particle::massFromPDGId(Particle::Proton)+0.135;
    const double m_proton = Particle::massFromPDGId( Particle::Proton ),
                 m2_proton = m_proton*m_proton,
                 m_min = m_proton+Particle::massFromPDGId( Particle::PiZero );

    const double mx = sqrt( mx2 );

    if ( mx < m_min || mx > 1.99 ) return StructureFunctions();

    int n_bin;
    double x_bin, dx;
    if ( mx < 1.11 ) {
      n_bin = 0;
      x_bin = mx-m_min;
      dx = 1.11-m_min; // Delta w bin sizes
    }
    else if ( mx < 1.77 ) { // w in [1.11, 1.77[
      dx = 0.015; // Delta w bin sizes
      n_bin = ( mx-1.11 )/dx + 1;
      x_bin = fmod( mx-1.11, dx );
    }
    else { // w in [1.77, 1.99[
      dx = 0.02; // Delta w bin sizes
      n_bin = ( mx-1.77 )/dx + 45;
      x_bin = fmod( mx-1.77, dx );
    }

    // values of a, b, c provided from the fits on ep data and retrieved from
    // http://dx.doi.org/10.1016/0550-3213(76)90231-5 with 1.110 <= w2 <=1.990
    const double abrass[56] = { 5.045, 5.126, 5.390,5.621, 5.913, 5.955,6.139,6.178,6.125, 5.999,
                                5.769, 5.622, 5.431,5.288, 5.175, 5.131,5.003,5.065,5.045, 5.078,
                                5.145, 5.156, 5.234,5.298, 5.371, 5.457,5.543,5.519,5.465, 5.384,
                                5.341, 5.320, 5.275,5.290, 5.330, 5.375,5.428,5.478,5.443, 5.390,
                                5.333, 5.296, 5.223,5.159, 5.146, 5.143,5.125,5.158,5.159, 5.178,
                                5.182, 5.195, 5.160,5.195, 5.163, 5.172 },
                 bbrass[56] = { 0.798, 1.052, 1.213,1.334,1.397,1.727,1.750,1.878,1.887,1.927,
                                2.041, 2.089, 2.148,2.205,2.344,2.324,2.535,2.464,2.564,2.610,
                                2.609, 2.678, 2.771,2.890,2.982,3.157,3.183,3.315,3.375,3.450,
                                3.477, 3.471, 3.554,3.633,3.695,3.804,3.900,4.047,4.290,4.519,
                                4.709, 4.757, 4.840,5.017,5.015,5.129,5.285,5.322,5.545,5.623,
                                5.775, 5.894, 6.138,6.151,6.301,6.542 },
                 cbrass[56] = { 0.043, 0.024, 0.000,-0.013,-0.023,-0.069,-0.060,-0.080,-0.065,-0.056,
                               -0.065,-0.056,-0.043,-0.034,-0.054,-0.018,-0.046,-0.015,-0.029,-0.048,
                               -0.032,-0.045,-0.084,-0.115,-0.105,-0.159,-0.164,-0.181,-0.203,-0.223,
                               -0.245,-0.254,-0.239,-0.302,-0.299,-0.318,-0.383,-0.393,-0.466,-0.588,
                               -0.622,-0.568,-0.574,-0.727,-0.665,-0.704,-0.856,-0.798,-1.048,-0.980,
                               -1.021,-1.092,-1.313,-1.341,-1.266,-1.473 };

    const double nu2 = pow( ( mx2-q2-m2_proton ) / ( 2.*m_proton ), 2 ),
                 logqq0 = log( ( nu2-q2 ) / pow( ( mx2-m2_proton ) / ( 2.*m_proton ), 2 ) ) / 2.,
                 gd2 = pow( 1. / ( 1-q2 / .71 ), 4 ); // dipole form factor of the proton

    const double sigLow = (n_bin == 0) ? 0. :
      exp( abrass[n_bin-1]+bbrass[n_bin-1]*logqq0+cbrass[n_bin-1]*pow( fabs( logqq0 ), 3 ) )*gd2;
    const double sigHigh =
      exp( abrass[n_bin]  +bbrass[n_bin]  *logqq0+cbrass[n_bin]  *pow( fabs( logqq0 ), 3 ) )*gd2;

    sigma_t = sigLow + x_bin*( sigHigh-sigLow )/dx;
    w1 = ( mx2-m2_proton )/( 8.*M_PI*M_PI*m_proton*Constants::alphaEM )/Constants::GeV2toBarn*1.e6 * sigma_t;
    w2 = w1 * q2/( q2-nu2 );

    return StructureFunctions( w1, w2 );
  }
}
