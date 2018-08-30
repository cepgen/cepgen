#include "CepGen/Physics/FormFactors.h"

#include "CepGen/Core/Exception.h"

#include "CepGen/Physics/ParticleProperties.h"
#include "CepGen/Physics/PDG.h"

#include "CepGen/StructureFunctions/FioreBrasse.h"
#include "CepGen/StructureFunctions/SuriYennie.h"
#include "CepGen/StructureFunctions/StructureFunctionsBuilder.h"

namespace CepGen
{
  const double FormFactors::mp_ = ParticleProperties::mass( PDG::proton );
  const double FormFactors::mp2_ = FormFactors::mp_*FormFactors::mp_;

  FormFactors
  FormFactors::trivial()
  {
    return FormFactors( 1.0, 1.0 );
  }

  FormFactors
  FormFactors::protonElastic( double q2 )
  {
    const double GE = pow( 1.+q2/0.71, -2. ), GE2 = GE*GE;
    const double GM = 2.79*GE, GM2 = GM*GM;
    return FormFactors( ( 4.*mp2_*GE2 + q2*GM2 ) / ( 4.*mp2_ + q2 ), GM2 );
  }

  FormFactors
  FormFactors::protonInelastic( double q2, double mi2, double mf2, StructureFunctions& sf )
  {
    switch ( sf.type ) {
      case SF::Type::ElasticProton:
        CG_WARNING( "FormFactors" ) << "Elastic proton form factors requested! Check your process definition!";
        return FormFactors::protonElastic( q2 );
      case SF::Type::SuriYennie:
        return FormFactors::suriYennie( q2, mi2, mf2 );
      case SF::Type::FioreBrasse:
        return FormFactors::fioreBrasse( q2, mi2, mf2 );
      default:
        return FormFactors::generic( q2, mi2, mf2, sf );
    }
  }

  FormFactors
  FormFactors::suriYennie( double q2, double mi2, double mf2 )
  {
    const double x = q2 / ( q2 + mf2 - mi2 );
    SF::SuriYennie suriyennie, sy = (SF::SuriYennie)suriyennie( x, q2 );
//std::cout << "---> " << sy.FM << "\t" << sy.F2*x/q2 << "\t" << sy.F2*x*sqrt(mi2)/q2 << std::endl;
    return FormFactors( sy.F2 * x * sqrt( mi2 ) / q2, sy.FM ); //FIXME
  }

  FormFactors
  FormFactors::fioreBrasse( double q2, double mi2, double mf2 )
  {
    const double x = q2 / ( q2 + mf2 - mi2 );
    SF::FioreBrasse fb, sf = fb( x, q2 );
    return FormFactors( sf.F2 * x / q2, -2.*sf.W1 / q2 );
  }

  FormFactors
  FormFactors::generic( double q2, double mi2, double mf2, StructureFunctions& sf )
  {
    const double x = q2 / ( q2 + mf2 - mi2 );
    sf = sf( x, q2 );
    sf.computeFL( x, q2 );
    return FormFactors( sf.F2 * x / q2, -2.*sf.F1( x, q2 ) / q2 );
  }

  double
  FormFactors::x( double q2, double w2, double m2 ) const
  {
    return 1./( 1.+( w2-mp2_ ) / q2+m2 );
  }

  std::ostream&
  operator<<( std::ostream& os, const FormFactors& ff )
  {
    os << Form( "Form factors: electric: Fe = %.3e ; magnetic: Fm = %.3e", ff.FE, ff.FM ).c_str();
    return os;
  }
}
