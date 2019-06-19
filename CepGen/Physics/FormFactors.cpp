#include "CepGen/Physics/FormFactors.h"
#include "CepGen/Physics/PDG.h"

#include "CepGen/Core/Exception.h"

#include "CepGen/StructureFunctions/SuriYennie.h"

namespace cepgen
{
  namespace ff
  {
    StandardDipole&
    StandardDipole::operator()( double q2 )
    {
      GE = pow( 1.+q2/0.71, -2. );
      GM = MU*GE;
      return *this;
    }
  }

  const double FormFactors::mp_ = PDG::get().mass( PDG::proton );
  const double FormFactors::mp2_ = FormFactors::mp_*FormFactors::mp_;

  FormFactors
  FormFactors::pointlikeScalar()
  {
    return FormFactors( 1.0, 0.0 );
  }

  FormFactors
  FormFactors::pointlikeFermion()
  {
    return FormFactors( 1.0, 1.0 );
  }

  FormFactors
  FormFactors::compositeScalar( double q2 )
  {
    throw CG_FATAL( "FormFactors:compositeScalar" ) << "Not yet supported!";
  }

  FormFactors
  FormFactors::protonElastic( double q2 )
  {
    const auto& fac = ff::StandardDipole()( q2 );
    const double GE2 = fac.GE*fac.GE, GM2 = fac.GM*fac.GM;
    return FormFactors( ( 4.*mp2_*GE2+q2*GM2 ) / ( 4.*mp2_+q2 ), GM2 );
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
        strfun::SuriYennie sy = strfun::SuriYennie()( xbj, q2 );
        return FormFactors( sy.F2 * xbj * sqrt( mi2 ) / q2, sy.FM ); //FIXME
      } break;
      default: {
        sf( xbj, q2 ).computeFL( xbj, q2 );
        return FormFactors( sf.F2 * xbj / q2, -2.*sf.F1( xbj, q2 ) / q2 );
      } break;
    }
  }

  std::ostream&
  operator<<( std::ostream& os, const FormFactors& ff )
  {
    return os << "FF{FE=" << ff.FE << ",FM=" << ff.FM << "}";
  }
}
