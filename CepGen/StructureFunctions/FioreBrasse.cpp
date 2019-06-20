#include "CepGen/StructureFunctions/FioreBrasse.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"

#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/Constants.h"

#include <complex>

namespace cepgen
{
  namespace strfun
  {
    FioreBrasse::Parameters
    FioreBrasse::Parameters::standard()
    {
      Parameters p;
      p.s0 = 1.14;
      p.norm = 0.021;
      p.resonances.emplace_back( Resonance{ -0.8377, 0.95, 0.1473, 1.0, 2.4617, 3./2. } ); // N*(1520)
      p.resonances.emplace_back( Resonance{ -0.37, 0.95, 0.1471, 0.5399, 2.4617, 5./2. } ); // N*(1680)
      p.resonances.emplace_back( Resonance{ 0.0038, 0.85, 0.1969, 4.2225, 1.5722, 3./2. } ); // Δ(1236)
      p.resonances.emplace_back( Resonance{ 0.5645, 0.1126, 1.3086, 19.2694, 4.5259, 1. } ); // exotic
      return p;
    }
    FioreBrasse::Parameters
    FioreBrasse::Parameters::alternative()
    {
      Parameters p;
      p.s0 = 1.2871;
      p.norm = 0.0207;
      p.resonances.emplace_back( Resonance{ -0.8070, 0.9632, 0.1387, 1.0, 2.6066, 3./2. } ); // N*(1520)
      p.resonances.emplace_back( Resonance{ -0.3640, 0.9531, 0.1239, 0.6086, 2.6066, 5./2. } ); // N*(1680)
      p.resonances.emplace_back( Resonance{ -0.0065, 0.8355, 0.2320, 4.7279, 1.4828, 3./2. } ); // Δ(1236)
      p.resonances.emplace_back( Resonance{ 0.5484, 0.1373, 1.3139, 14.7267, 4.6041, 1. } ); // exotic
      return p;
    }

    FioreBrasse::FioreBrasse( const ParametersList& params ) :
      Parameterisation( params )
    {
      const auto& model = params.get<std::string>( "model", "standard" );
      if ( model == "standard" )
        params_ = Parameters::standard();
      else if ( model == "alternative" )
        params_ = Parameters::alternative();
      else
        throw CG_FATAL( "FioreBrasse" ) << "Invalid modelling selected: " << model << "!";
    }

    FioreBrasse&
    FioreBrasse::operator()( double xbj, double q2 )
    {
      std::pair<double,double> nv = { xbj, q2 };
      if ( nv == old_vals_ )
        return *this;
      old_vals_ = nv;

      const double akin = 1. + 4.*mp2_ * xbj*xbj/q2;
      const double prefactor = q2*( 1.-xbj ) / ( 4.*M_PI*constants::ALPHA_EM*akin );
      const double s = q2*( 1.-xbj )/xbj + mp2_;

      double ampli_res = 0., ampli_bg = 0., ampli_tot = 0.;
      for ( unsigned short i = 0; i < 3; ++i ) { //FIXME 4??
        const Parameters::Resonance& res = params_.resonances[i];
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
        const Parameters::Resonance& res = params_.resonances[3];
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

      F2 = prefactor*ampli_tot;
      return *this;
    }
  }
}

REGISTER_STRFUN( FioreBrasse, strfun::FioreBrasse )
