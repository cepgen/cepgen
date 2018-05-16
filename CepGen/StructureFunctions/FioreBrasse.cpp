#include "CepGen/StructureFunctions/FioreBrasse.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/ParticleProperties.h"
#include "CepGen/Physics/Constants.h"

#include <complex>

namespace CepGen
{
  namespace SF
  {
    FioreBrasse::Parameterisation
    FioreBrasse::Parameterisation::standard()
    {
      Parameterisation p;
      p.s0 = 1.14;
      p.norm = 0.021;
      p.resonances.emplace_back( -0.8377, 0.95, 0.1473, 1.0, 2.4617, 3./2. );
      p.resonances.emplace_back( -0.37, 0.95, 0.1471, 0.5399, 2.4617, 5./2. );
      p.resonances.emplace_back( 0.0038, 0.85, 0.1969, 4.2225, 1.5722, 3./2. );
      p.resonances.emplace_back( 0.5645, 0.1126, 1.3086, 19.2694, 4.5259, 1. );
      return p;
    }
    FioreBrasse::Parameterisation
    FioreBrasse::Parameterisation::alternative()
    {
      Parameterisation p;
      p.s0 = 1.2871;
      p.norm = 0.0207;
      p.resonances.emplace_back( -0.8070, 0.9632, 0.1387, 1.0, 2.6066, 3./2. );
      p.resonances.emplace_back( -0.3640, 0.9531, 0.1239, 0.6086, 2.6066, 5./2. );
      p.resonances.emplace_back( -0.0065, 0.8355, 0.2320, 4.7279, 1.4828, 3./2. );
      p.resonances.emplace_back( 0.5484, 0.1373, 1.3139, 14.7267, 4.6041, 1. );
      return p;
    }

    FioreBrasse::FioreBrasse( const FioreBrasse::Parameterisation& params ) :
      W1( 0. ), W2( 0. ), params_( params )
    {}

    FioreBrasse
    FioreBrasse::operator()( double q2, double xbj ) const
    {
      const double akin = 1. + 4.*mp2_ * xbj*xbj/q2;
      const double prefactor = q2*( 1.-xbj ) / ( 4.*M_PI*Constants::alphaEM*akin );
      const double s = q2*( 1.-xbj )/xbj + mp2_;

      double ampli_res = 0., ampli_bg = 0., ampli_tot = 0.;
      for ( unsigned short i = 0; i < 3; ++i ) { //FIXME 4??
        const Parameterisation::ResonanceParameters res = params_.resonances[i];
        if ( !res.enabled )
          continue;
        const double sqrts0 = sqrt( params_.s0 );

        std::complex<double> alpha;
        if ( s > params_.s0 )
          alpha = std::complex<double>( res.alpha0 + res.alpha2*sqrts0 + res.alpha1*s, res.alpha2*sqrt( s-params_.s0 ) );
        else
          alpha = std::complex<double>( res.alpha0 + res.alpha1*s + res.alpha2*( sqrts0 - sqrt( params_.s0 - s ) ), 0. );

        double formfactor = 1./pow( 1. + q2/res.q02, 2 );
        double denom = pow( res.spin-std::real( alpha ), 2 ) + pow( std::imag( alpha ), 2 );
        double ampli_imag = res.a*formfactor*formfactor*std::imag( alpha )/denom;
        ampli_res += ampli_imag;
      }
      {
        const Parameterisation::ResonanceParameters res = params_.resonances[3];
        double sE = res.alpha2, sqrtsE = sqrt( sE );
        std::complex<double> alpha;
        if ( s > sE )
          alpha = std::complex<double>( res.alpha0 + res.alpha1*sqrtsE, res.alpha1*sqrt( s-sE ) );
        else
          alpha = std::complex<double>( res.alpha0 + res.alpha1*( sqrtsE - sqrt( sE-s ) ), 0. );
        double formfactor = 1./pow( 1. + q2/res.q02, 2 );
        double sp = 1.5*res.spin;
        double denom = pow( sp-std::real( alpha ), 2 ) + pow( std::imag( alpha ), 2 );
        ampli_bg = res.a*formfactor*formfactor*std::imag( alpha )/denom;
      }
      ampli_tot = params_.norm*( ampli_res+ampli_bg );

      CG_DEBUG_LOOP( "FioreBrasse:amplitudes" )
        << "Amplitudes:\n\t"
        << " resonance part:  " << ampli_res << ",\n\t"
        << " background part: " << ampli_bg << ",\n\t"
        << " total (with norm.): " << ampli_tot << ".";

      FioreBrasse fb;
      fb.F2 = prefactor*ampli_tot;
      fb.computeFL( q2, xbj );
      return fb;
    }

    FioreBrasse
    FioreBrasse::operator()( double q2, double xbj, bool ) const
    {
      const double m_min = mp_+ParticleProperties::mass( PDG::PiZero );

      const double mx2 = mp2_ + q2*( 1.-xbj )/xbj, mx = sqrt( mx2 );

      FioreBrasse fb;
      if ( mx < m_min || mx > 1.99 ) {
        CG_WARNING( "FioreBrasse" )
          << "Fiore-Brasse form factors to be retrieved for an invalid MX value:\n\t"
          << mx << " GeV, while allowed range is [1.07, 1.99] GeV.";
        return fb;
      }

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
      const double a[56] = { 5.045, 5.126, 5.390,5.621, 5.913, 5.955,6.139,6.178,6.125, 5.999,
                             5.769, 5.622, 5.431,5.288, 5.175, 5.131,5.003,5.065,5.045, 5.078,
                             5.145, 5.156, 5.234,5.298, 5.371, 5.457,5.543,5.519,5.465, 5.384,
                             5.341, 5.320, 5.275,5.290, 5.330, 5.375,5.428,5.478,5.443, 5.390,
                             5.333, 5.296, 5.223,5.159, 5.146, 5.143,5.125,5.158,5.159, 5.178,
                             5.182, 5.195, 5.160,5.195, 5.163, 5.172 },
                   b[56] = { 0.798, 1.052, 1.213,1.334,1.397,1.727,1.750,1.878,1.887,1.927,
                             2.041, 2.089, 2.148,2.205,2.344,2.324,2.535,2.464,2.564,2.610,
                             2.609, 2.678, 2.771,2.890,2.982,3.157,3.183,3.315,3.375,3.450,
                             3.477, 3.471, 3.554,3.633,3.695,3.804,3.900,4.047,4.290,4.519,
                             4.709, 4.757, 4.840,5.017,5.015,5.129,5.285,5.322,5.545,5.623,
                             5.775, 5.894, 6.138,6.151,6.301,6.542 },
                   c[56] = { 0.043, 0.024, 0.000,-0.013,-0.023,-0.069,-0.060,-0.080,-0.065,-0.056,
                            -0.065,-0.056,-0.043,-0.034,-0.054,-0.018,-0.046,-0.015,-0.029,-0.048,
                            -0.032,-0.045,-0.084,-0.115,-0.105,-0.159,-0.164,-0.181,-0.203,-0.223,
                            -0.245,-0.254,-0.239,-0.302,-0.299,-0.318,-0.383,-0.393,-0.466,-0.588,
                            -0.622,-0.568,-0.574,-0.727,-0.665,-0.704,-0.856,-0.798,-1.048,-0.980,
                            -1.021,-1.092,-1.313,-1.341,-1.266,-1.473 };

      const double d = 3.0;
      const double nu = 0.5 * ( q2 + mx2 - mp2_ ) / mp_, nu2 = nu*nu,
                   logqq0 = 0.5 * log( ( nu2+q2 ) / pow( ( mx2-mp2_ ) / ( 2.*mp_ ), 2 ) );
      const double gd2 = pow( 1. / ( 1+q2 / .71 ), 4 ); // dipole form factor of the proton

      const double sigLow = ( n_bin == 0 ) ? 0. :
        gd2 * exp( a[n_bin-1] + b[n_bin-1]*logqq0 + c[n_bin-1]*pow( fabs( logqq0 ), d ) );
      const double sigHigh =
        gd2 * exp( a[n_bin]   + b[n_bin]  *logqq0 + c[n_bin]  *pow( fabs( logqq0 ), d ) );

      const double sigma_t = sigLow + x_bin*( sigHigh-sigLow )/dx;
      const double w1 = ( mx2-mp2_ )/( 8.*M_PI*M_PI*mp_*Constants::alphaEM )/Constants::GeV2toBarn*1.e6 * sigma_t;
      const double w2 = w1 * q2 / ( q2+nu2 );

      fb.W1 = w1; //FIXME
      fb.W2 = w2; //FIXME
      return fb;
    }
  }
}
