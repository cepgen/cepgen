#include "SuriYennie.h"

namespace CepGen
{
  namespace SF
  {
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
      sy.F1 = 0.5*sy.F2/xbj; // Callan-Gross relation FIXME
      return sy;
    }
  }
}
