#include "SuriYennie.h"

namespace CepGen
{
  namespace SF
  {
    SuriYennieParameterisation
    SuriYennieParameterisation::standard()
    {
      SuriYennieParameterisation p;
      p.C1 = 0.86926;
      p.C2 = 2.23422;
      p.D1 = 0.12549;
      p.rho2 = 0.585;
      p.Cp = 0.96;
      p.Bp = 0.63;
      return p;
    }

    SuriYennieParameterisation
    SuriYennieParameterisation::alternative()
    {
      SuriYennieParameterisation p;
      p.C1 = 0.6303;
      p.C2 = 2.3049;
      p.D1 = 0.04681;
      p.rho2 = 1.05;
      p.Cp = 1.23;
      p.Bp = 0.61;
      return p;
    }

    StructureFunctions
    SuriYennie( double q2, double xbj, const SuriYennieParameterisation& param )
    {
      const double mp = Particle::massFromPDGId( Particle::Proton ), mp2 = mp*mp;
      const double mx2 = q2 * ( 1.-xbj )/xbj + mp2, // [GeV^2]
                   nu = 0.5 * ( q2 + mx2 - mp2 ) / mp; // [GeV]
      const double dm2 = mx2-mp2, Xpr = q2/( q2+mx2 ), En = dm2+q2, Tau = 0.25 * q2/mp2, MQ = param.rho2+q2;

      StructureFunctions sy;
      sy.FM = ( param.C1*dm2*pow( param.rho2/MQ, 2 ) + ( param.C2*mp2*pow( 1.-Xpr, 4 ) ) / ( 1.+Xpr*( Xpr*param.Cp-2.*param.Bp ) ) )/q2;
      const double FE = ( Tau*sy.FM + param.D1*dm2*q2*param.rho2/mp2*pow( dm2/MQ/En, 2 ) ) / ( 1.+0.25*En*En/mp2/q2 );

      const double w2 = 2.*mp*FE, w1 = 0.5 * sy.FM*q2/mp;

      sy.F2 = nu/mp*w2;
      return sy;
    }
  }
}
