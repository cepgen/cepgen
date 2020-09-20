#include "CepGen/StructureFunctions/SuriYennie.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"

#include "CepGen/Core/Exception.h"

#include <math.h>

namespace cepgen
{
  namespace strfun
  {
    SuriYennie::Parameters
    SuriYennie::Parameters::standard()
    {
      Parameters p;
      p.C1 = 0.86926;
      p.C2 = 2.23422;
      p.D1 = 0.12549;
      p.rho2 = 0.585;
      p.Cp = 0.96;
      p.Bp = 0.63;
      return p;
    }

    SuriYennie::Parameters
    SuriYennie::Parameters::alternative()
    {
      Parameters p;
      p.C1 = 0.6303;
      p.C2 = 2.3049;
      p.D1 = 0.04681;
      p.rho2 = 1.05;
      p.Cp = 1.23;
      p.Bp = 0.61;
      return p;
    }

    SuriYennie::SuriYennie( const ParametersList& params ) :
      Parameterisation( params ),
      W1( 0. ), W2( 0. ), FE( 0. ), FM( 0. )
    {
      const auto& model = params.get<std::string>( "model", "standard" );
      if ( model == "standard" )
        params_ = Parameters::standard();
      else if ( model == "alternative" )
        params_ = Parameters::alternative();
      else
        throw CG_FATAL( "SuriYennie" ) << "Invalid modelling selected: " << model << "!";
    }

    SuriYennie&
    SuriYennie::operator()( double xbj, double q2 )
    {
      std::pair<double,double> nv = { xbj, q2 };
      if ( nv == old_vals_ )
        return *this;
      old_vals_ = nv;

      const double inv_q2 = 1./q2;
      const double dm2 = q2*( 1.-xbj )/xbj;
      const double mx2 = mp2_ + dm2; // [GeV^2]
      const double en = dm2+q2, en2 = en*en; // [GeV^2]
      const double nu = 0.5 * en / mp_, Xpr = q2/( q2+mx2 ), tau = 0.25 * q2/mp2_;
      const double mq = params_.rho2+q2;

      FM = ( params_.C1*dm2*pow( params_.rho2/mq, 2 )
            + ( params_.C2*mp2_*pow( 1.-Xpr, 4 ) ) / ( 1.+Xpr*( Xpr*params_.Cp-2.*params_.Bp ) ) ) * inv_q2;
      FE = ( tau*FM + params_.D1*dm2*q2*params_.rho2/mp2_*pow( dm2/mq, 2 )/en2 ) / ( 1.+0.25*en2/mp2_*inv_q2 );

      //const double w2 = 2.*mp*FE;

      W1 = 0.5*FM*q2/mp_;
      W2 = 2.*mp_*FE;
      F2 = 2.*nu*FE;
      return *this;
    }
  }
}

REGISTER_STRFUN( SuriYennie, strfun::SuriYennie )
