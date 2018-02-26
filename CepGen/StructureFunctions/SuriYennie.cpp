#include "SuriYennie.h"
#include "CepGen/Physics/ParticleProperties.h"
#include <math.h>

namespace CepGen
{
  namespace SF
  {
    SuriYennie::Parameterisation
    SuriYennie::Parameterisation::standard()
    {
      Parameterisation p;
      p.C1 = 0.86926;
      p.C2 = 2.23422;
      p.D1 = 0.12549;
      p.rho2 = 0.585;
      p.Cp = 0.96;
      p.Bp = 0.63;
      return p;
    }

    SuriYennie::Parameterisation
    SuriYennie::Parameterisation::alternative()
    {
      Parameterisation p;
      p.C1 = 0.6303;
      p.C2 = 2.3049;
      p.D1 = 0.04681;
      p.rho2 = 1.05;
      p.Cp = 1.23;
      p.Bp = 0.61;
      return p;
    }

    SuriYennie::SuriYennie( const SuriYennie::Parameterisation& param ) :
      F1( 0. ), FE( 0. ), FM( 0. ), params_( param )
    {}

    SuriYennie
    SuriYennie::operator()( double q2, double xbj ) const
    {
      const double mp = ParticleProperties::mass( Proton ), mp2 = mp*mp;
      const double inv_q2 = 1./q2;
      const double dm2 = q2*( 1.-xbj )/xbj;
      const double mx2 = mp2 + dm2; // [GeV^2]
      const double en = dm2+q2, en2 = en*en; // [GeV^2]
      const double nu = 0.5 * en / mp, Xpr = q2/( q2+mx2 ), tau = 0.25 * q2/mp2;
      const double mq = params_.rho2+q2;

      SuriYennie sy;
      sy.FM = ( params_.C1*dm2*pow( params_.rho2/mq, 2 )
            + ( params_.C2*mp2*pow( 1.-Xpr, 4 ) ) / ( 1.+Xpr*( Xpr*params_.Cp-2.*params_.Bp ) ) ) * inv_q2;
      sy.FE = ( tau*sy.FM + params_.D1*dm2*q2*params_.rho2/mp2*pow( dm2/mq, 2 )/en2 ) / ( 1.+0.25*en2/mp2*inv_q2 );

      //const double w2 = 2.*mp*sy.FE;

      sy.F1 = 0.5*sy.FM*q2/mp;
      sy.F2 = 2.*nu*sy.FE;
      return sy;
    }
  }
}
