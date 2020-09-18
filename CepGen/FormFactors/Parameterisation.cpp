#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/StructureFunctions/SuriYennie.h"

#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/PDG.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/String.h"

#include <cmath>
#include <cassert>

namespace cepgen
{
  namespace formfac
  {
    Parameterisation::Parameterisation() :
      NamedModule<std::string>( ParametersList() ),
      mp_( PDG::get().mass( PDG::proton ) ), mp2_( mp_*mp_ ),
      last_q2_( -1. ),
      FE( 0. ), FM( 0. ), GE( 0. ), GM( 0. )
    {}

    Parameterisation::Parameterisation( const Parameterisation& param ) :
      NamedModule<std::string>( param.parameters() ),
      mp_( param.mp_ ), mp2_( param.mp2_ ),
      last_q2_( -1. ),
      FE( param.FE ), FM( param.FM ), GE( param.GE ), GM( param.GM )
    {}

    Parameterisation::Parameterisation( const ParametersList& params ) :
      NamedModule<std::string>( params ),
      mp_( PDG::get().mass( PDG::proton ) ), mp2_( mp_*mp_ ),
      last_q2_( -1. ),
      FE( 0. ), FM( 0. ), GE( 0. ), GM( 0. )
    {}

    double
    Parameterisation::tau( double q2 ) const
    {
      if ( mp2_ <= 0. )
        throw CG_FATAL( "FormFactors:tau" )
          << "Invalid proton mass! check the form factors constructor!";
      return 0.25*q2/mp2_;
    }

    Parameterisation&
    Parameterisation::operator()( const mode::Beam& type, double q2, double mi2, double mf2 )
    {
      last_q2_ = q2;
      switch ( type ) {
        case mode::Beam::invalid:
        case mode::Beam::CompositeScalar:
          throw CG_FATAL( "FormFactors" )
            << type << " mode is not yet supported!";
        case mode::Beam::PointLikeScalar:
          FE = 1., FM = 0.; break;
        case mode::Beam::PointLikeFermion:
          FE = FM = 1.; break;
        case mode::Beam::ProtonElastic: {
          compute( q2 );
          const double GE2 = GE*GE, GM2 = GM*GM;
          FE = ( 4.*mp2_*GE2+q2*GM2 ) / ( 4.*mp2_+q2 );
          FM = GM2;
        } break;
        case mode::Beam::ProtonInelastic: {
          const double xbj = q2/( q2+mf2-mi2 );
          if ( !str_fun_ )
            throw CG_FATAL( "FormFactors" )
              << "Inelastic proton form factors computation requires "
              << "a structure functions definition!";
          switch ( (strfun::Type)str_fun_->name() ) {
            case strfun::Type::ElasticProton:
              throw CG_FATAL( "FormFactors" )
                << "Elastic proton form factors requested!\n"
                << "Check your process definition!";
            case strfun::Type::SuriYennie: { //FIXME
              static strfun::SuriYennie sy;
              sy = sy( xbj, q2 );
              FE = sy.F2 * xbj * sqrt( mi2 ) / q2;
              FM = sy.FM;
            } break;
            default: {
              ( *str_fun_ )( xbj, q2 ).computeFL( xbj, q2 );
              FE = str_fun_->F2 * xbj / q2;
              FM = -2.*str_fun_->F1( xbj, q2 ) / q2;
            } break;
          }
        } break;
      }
      return *this;
    }

    //------------------------------------------------------------------

    std::ostream&
    operator<<( std::ostream& os, const Parameterisation* ff )
    {
      os << ff->name();
      if ( ff->last_q2_ >= 0. )
        os << "(Q²=" << ff->last_q2_ << " GeV²): "
           << "FE=" << ff->FE << ",FM=" << ff->FM;
      return os;
    }

    std::ostream&
    operator<<( std::ostream& os, const Parameterisation& ff )
    {
      return os << &ff;
    }
  }
}
