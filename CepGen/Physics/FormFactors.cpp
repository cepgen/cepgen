#include "CepGen/Physics/FormFactors.h"

#include "CepGen/Core/Exception.h"

#include "CepGen/Physics/ParticleProperties.h"
#include "CepGen/Physics/PDG.h"

#include "CepGen/StructureFunctions/SuriYennie.h"

namespace cepgen
{
  const double FormFactors::mp_ = particleproperties::mass( PDG::proton );
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
  FormFactors::protonInelastic( double q2, double mi2, double mf2, strfun::Parameterisation& sf )
  {
    const double xbj = q2 / ( q2 + mf2 - mi2 );
    switch ( sf.type ) {
      case strfun::Type::ElasticProton:
        CG_WARNING( "FormFactors" ) << "Elastic proton form factors requested! Check your process definition!";
        return FormFactors::protonElastic( q2 );
      case strfun::Type::SuriYennie: {
        strfun::SuriYennie suriyennie, sy = (strfun::SuriYennie)suriyennie( xbj, q2 );
        return FormFactors( sy.F2 * xbj * sqrt( mi2 ) / q2, sy.FM ); //FIXME
      } break;
      default: {
        sf = sf( xbj, q2 );
        sf.computeFL( xbj, q2 );
        return FormFactors( sf.F2 * xbj / q2, -2.*sf.F1( xbj, q2 ) / q2 );
      } break;
    }
  }

  std::ostream&
  operator<<( std::ostream& os, const FormFactors& ff )
  {
    os << Form( "Form factors: electric: Fe = %.3e ; magnetic: Fm = %.3e", ff.FE, ff.FM ).c_str();
    return os;
  }
}
