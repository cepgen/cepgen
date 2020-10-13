#include "CepGen/StructureFunctions/SuriYennie.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/Physics.h"

#include <math.h>

namespace cepgen
{
  namespace strfun
  {
    SuriYennie::SuriYennie( const ParametersList& params ) :
      Parameterisation( params ),
      W1( 0. ), W2( 0. ), FE( 0. ), FM( 0. )
    {
      const auto& model = params.get<std::string>( "model", "standard" );
      if ( model == "standard" ) {
        params_.C1 = 0.86926;
        params_.C2 = 2.23422;
        params_.D1 = 0.12549;
        params_.rho2 = 0.585;
        params_.Cp = 0.96;
        params_.Bp = 0.63;
      }
      else if ( model == "alternative" ) {
        params_.C1 = 0.6303;
        params_.C2 = 2.3049;
        params_.D1 = 0.04681;
        params_.rho2 = 1.05;
        params_.Cp = 1.23;
        params_.Bp = 0.61;
      }
      else { // custom model
        params_.C1 = params.get<double>( "C1" );
        params_.C2 = params.get<double>( "C2" );
        params_.D1 = params.get<double>( "D1" );
        params_.rho2 = params.get<double>( "rho2" );
        params_.Cp = params.get<double>( "Cp" );
        params_.Bp = params.get<double>( "Bp" );
      }
    }

    SuriYennie&
    SuriYennie::eval( double xbj, double q2 )
    {
      const double mx2 = utils::mX2( xbj, q2, mp2_ ), dm2 = mx2-mp2_; // [GeV^2]
      const double en = q2+dm2, en2 = en*en; // [GeV^2]
      const double nu = 0.5 * en / mp_, x_pr = q2/( q2+mx2 ), tau = 0.25 * q2/mp2_;
      const double mq = params_.rho2+q2;

      const double inv_q2 = 1./q2;

      FM = ( params_.C1*dm2*pow( params_.rho2/mq, 2 )
            + ( params_.C2*mp2_*pow( 1.-x_pr, 4 ) ) / ( 1.+x_pr*( x_pr*params_.Cp-2.*params_.Bp ) ) ) * inv_q2;
      FE = ( tau*FM + params_.D1*dm2*q2*params_.rho2/mp2_*pow( dm2/mq, 2 )/en2 ) / ( 1.+0.25*en2/mp2_*inv_q2 );

      W1 = 0.5*FM*q2/mp_;
      W2 = 2.*mp_*FE;
      F2 = 2.*nu*FE;
      return *this;
    }
  }
}

REGISTER_STRFUN( SuriYennie, strfun::SuriYennie )
