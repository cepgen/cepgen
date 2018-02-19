#include "FormFactors.h"

#include "CepGen/Core/Exception.h"

#include "CepGen/Physics/ParticleProperties.h"

#include "CepGen/StructureFunctions/ALLM.h"
#include "CepGen/StructureFunctions/BlockDurandHa.h"
#include "CepGen/StructureFunctions/FioreBrasse.h"
#include "CepGen/StructureFunctions/GenericLHAPDF.h"
#include "CepGen/StructureFunctions/SuriYennie.h"
#include "CepGen/StructureFunctions/SzczurekUleshchenko.h"

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
    const double mp2 = ParticleProperties::mass( Proton )*ParticleProperties::mass( Proton );
    const double GE = pow( 1.+q2/0.71, -2. ), GE2 = GE*GE;
    const double GM = 2.79*GE, GM2 = GM*GM;
    return FormFactors( ( 4.*mp2*GE2 + q2*GM2 ) / ( 4.*mp2 + q2 ), GM2 );
  }

  FormFactors
  FormFactors::ProtonInelastic( const StructureFunctions::Type& sf, double q2, double mi2, double mf2 )
  {
    switch ( sf ) {
      case StructureFunctions::ElasticProton:
        InWarning( "Elastic proton form factors requested! Check your process definition!" );
        return FormFactors::ProtonElastic( q2 );
      case StructureFunctions::SuriYennie:
        return FormFactors::SuriYennie( q2, mi2, mf2 );
      case StructureFunctions::SzczurekUleshchenko:
        return FormFactors::SzczurekUleshchenko( q2, mi2, mf2 );
      case StructureFunctions::FioreBrasse:
        return FormFactors::FioreBrasse( q2, mi2, mf2 );
      default: throw Exception( __PRETTY_FUNCTION__, "Invalid structure functions required!", FatalError );
    }
  }

  FormFactors
  FormFactors::SuriYennie( double q2, double mi2, double mf2 )
  {
    const double x = q2 / ( q2 + mf2 - mi2 );
    const SF::SuriYennie suriyennie, sy = suriyennie( q2, x );
//std::cout << "---> " << sy.FM << "\t" << sy.F2*x/q2 << "\t" << sy.F2*x*sqrt(mi2)/q2 << std::endl;
    return FormFactors( sy.F2 * x * sqrt( mi2 ) / q2, sy.FM ); //FIXME
  }

  FormFactors
  FormFactors::FioreBrasse( double q2, double mi2, double mf2 )
  {
    const double x = q2 / ( q2 + mf2 - mi2 );
    SF::FioreBrasse fb, sf = fb( q2, x );
    return FormFactors( sf.F2 * x / q2, -2.*sf.W1 / q2 );
  }

  FormFactors
  FormFactors::SzczurekUleshchenko( double q2, double mi2, double mf2 )
  {
    const double x = q2 / ( q2 + mf2 - mi2 );
    SF::SzczurekUleshchenko su, sf = su( q2, x );
    return FormFactors( sf.F2 * x / q2, -2.*sf.F1 / q2 );
  }

  double
  FormFactors::x( double q2, double w2, double m2 ) const
  {
    const double mp2 = ParticleProperties::mass( Proton )*ParticleProperties::mass( Proton );
    return 1./( 1.+( w2-mp2 ) / q2+m2 );
  }

  std::ostream&
  operator<<( std::ostream& os, const FormFactors& ff )
  {
    os << Form( "Form factors: electric: Fe = %.3e ; magnetic: Fm = %.3e", ff.FE, ff.FM ).c_str();
    return os;
  }
}
