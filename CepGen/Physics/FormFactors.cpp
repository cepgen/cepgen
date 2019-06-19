#include "CepGen/Physics/FormFactors.h"
#include "CepGen/Physics/PDG.h"

#include "CepGen/Core/Exception.h"

#include "CepGen/StructureFunctions/SuriYennie.h"

namespace cepgen
{
  namespace ff
  {
    const double Parameterisation::mp_ = PDG::get().mass( PDG::proton );
    const double Parameterisation::mp2_ = Parameterisation::mp_*Parameterisation::mp_;

    std::shared_ptr<Parameterisation>
    Parameterisation::build( const ParametersList& params )
    {
      const auto model = (Model)params.get<int>( "model" );
      switch ( model ) {
        case Model::StandardDipole: return std::make_shared<StandardDipole>( params );
        //case Model::ArringtonEtAl: return std::make_shared<ArringtonEtAl>( params );
        //case Model::BrashEtAl: return std::make_shared<BrashEtAl>( params );
        default:
          throw CG_FATAL( "FormFactors" ) << "Invalid FF modelling requested: " << (int)model << "!";
      }
    }

    Parameterisation::Parameterisation( const ParametersList& params ) :
      model_( (Model)params.get<int>( "model", (int)Model::Invalid ) ),
      type_( (Type)params.get<int>( "type", (int)Type::Invalid ) ),
      str_fun_( strfun::Parameterisation::build( params.get<ParametersList>( "structureFunctions" ) ) )
    {}

    Parameterisation&
    Parameterisation::operator()( double q2, double mi2, double mf2 )
    {
      switch ( type_ ) {
        case Type::Invalid:
        case Type::CompositeScalar:
          throw CG_FATAL( "FormFactors" ) << "Not yet supported!";
        case Type::PointLikeScalar:
          FE = 1.;
          FM = 0.;
          break;
        case Type::PointLikeFermion:
          FE = 1.;
          FM = 1.;
          break;
        case Type::ProtonElastic: {
          compute( q2 );
          const double GE2 = GE*GE, GM2 = GM*GM;
          FE = ( 4.*mp2_*GE2+q2*GM2 ) / ( 4.*mp2_+q2 );
          FM = GM2;
        } break;
        case Type::ProtonInelastic: {
          const double xbj = q2 / ( q2 + mf2 - mi2 );
          switch ( str_fun_->type ) {
            case strfun::Type::ElasticProton:
              throw CG_FATAL( "FormFactors" )
                << "Elastic proton form factors requested! Check your process definition!";
            case strfun::Type::SuriYennie: { //FIXME
              strfun::SuriYennie sy = strfun::SuriYennie()( xbj, q2 );
              FE = sy.F2 * xbj * sqrt( mi2 ) / q2;
              FM = sy.FM;
            } break;
            default: {
              str_fun_->operator()( xbj, q2 ).computeFL( xbj, q2 );
              FE = str_fun_->F2 * xbj / q2;
              FM = -2.*str_fun_->F1( xbj, q2 ) / q2;
            } break;
          }
        } break;
      }
      return *this;
    }

    void
    StandardDipole::compute( double q2 )
    {
      GE = pow( 1.+q2/0.71, -2. );
      GM = MU*GE;
    }
  }

  std::ostream&
  operator<<( std::ostream& os, const ff::Parameterisation& formfac )
  {
    return os << "FF{FE=" << formfac.FE << ",FM=" << formfac.FM << "}";
  }
}
