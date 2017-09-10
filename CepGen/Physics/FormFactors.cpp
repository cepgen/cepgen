#include "FormFactors.h"

namespace CepGen
{
  FormFactors
  FormFactors::Trivial()
  {
    return FormFactors( 1.0, 1.0 );
  }

  FormFactors
  FormFactors::ProtonElastic( double q2 )
  {
    const double mp = Particle::massFromPDGId( Particle::Proton ), mp2 = mp*mp;
    const double GE = pow( 1.+q2/0.71, -2. ), GE2 = GE*GE;
    const double GM = 2.79*GE, GM2 = GM*GM;
    return FormFactors( ( 4.*mp2*GE2 + q2*GM2 ) / ( 4.*mp2 + q2 ), GM2 );
  }

  FormFactors
  FormFactors::ProtonInelastic( const StructureFunctionsType& sf, double q2, double mi2, double mf2 )
  {
    switch ( sf ) {
      case StructureFunctionsType::ElasticProton:
        InWarning( "Elastic proton form factors requested! Check your process definition!" );
        return FormFactors::ProtonElastic( q2 );
      case StructureFunctionsType::SuriYennie:
        return FormFactors::SuriYennie( q2, mi2, mf2 );
      case StructureFunctionsType::SzczurekUleshchenko:
        return FormFactors::SzczurekUleshchenko( q2, mi2, mf2 );
      case StructureFunctionsType::Fiore:
      case StructureFunctionsType::FioreSea:
      case StructureFunctionsType::FioreVal:
        return FormFactors::FioreBrasse( q2, mi2, mf2 );
      default: throw Exception( __PRETTY_FUNCTION__, "Invalid structure functions required!", FatalError );
    }
  }

  FormFactors
  FormFactors::SuriYennie( double q2, double mi2, double mf2 )
  {
    struct SuriYennieParameters {
      SuriYennieParameters( double rho, double bp, double cp, double cc1, double cc2, double dd1 ) :
        rho( rho ), Bp( bp ), Cp( cp ), cc1( cc1 ), cc2( cc2 ), dd1( dd1 ) {}
      double rho; // in GeV**2
      double Bp, Cp;
      double cc1, cc2, dd1;
    };
    // values extracted from experimental fits
    SuriYennieParameters sy( 0.585, 0.63, 0.96, 0.86926, 2.23422, 0.12549 );
    //SuriYennieParameters sy( 1.05, 0.61, 1.23, 0.6303, 2.2049, 0.0468 );
    const double x = q2 / ( q2 + mf2 ),
                 dm2 = mf2-mi2,
                 en = q2 + dm2,
                 tau = q2 / mi2 / 4.,
                 rhot = sy.rho+q2;
    const double rho_norm = sy.rho/rhot;

    FormFactors ff;
    ff.FM = ( -1./q2 ) * ( -sy.cc1*rho_norm*rho_norm*dm2 - sy.cc2*mi2*pow( 1.-x, 4 )/( x*( x*sy.Cp-2*sy.Bp )+1. ) );
    ff.FE = ( tau*ff.FM + sy.dd1*dm2*q2*rho_norm*pow( dm2/en, 2 )/( rhot*mi2 ) )/( 1. + en*en/( 4.*mi2*q2 ) );
    return ff;
  }

  FormFactors
  FormFactors::FioreBrasse( double q2, double mi2, double mf2 )
  {
    const double x = q2 / ( q2 + mf2 - mi2 ), k = 2.*sqrt( mi2 );
    StructureFunctions sf = StructureFunctions::FioreBrasse( q2, x );
    return FormFactors( sf.F2 / k, -sf.F1*k / q2 );
  }

  FormFactors
  FormFactors::SzczurekUleshchenko( double q2, double mi2, double mf2 )
  {
    const double x = q2 / ( q2 + mf2 - mi2 );
    StructureFunctions sf = StructureFunctions::SzczurekUleshchenko( q2, x );
    return FormFactors( sf.F2 * x / q2, -2.*sf.F1 / q2 );
  }

  std::ostream&
  operator<<( std::ostream& os, const FormFactors& ff )
  {
    os << Form( "Form factors: electric: Fe = %.3e ; magnetic: Fm = %.3e", ff.FE, ff.FM ).c_str();
    return os;
  }
}

