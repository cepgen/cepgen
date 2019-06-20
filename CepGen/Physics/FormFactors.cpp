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
        case Model::ArringtonEtAl: return std::make_shared<ArringtonEtAl>( params );
        case Model::BrashEtAl: return std::make_shared<BrashEtAl>( params );
        default:
          throw CG_FATAL( "FormFactors" ) << "Invalid FF modelling requested: " << (int)model << "!";
      }
    }

    Parameterisation::Parameterisation() :
      model_( Model::Invalid ), type_( Type::Invalid ),
      FE( 0. ), FM( 0. ), GE( 0. ), GM( 0. )
    {}

    Parameterisation::Parameterisation( const ParametersList& params ) :
      model_( (Model)params.get<int>( "model", (int)Model::Invalid ) ),
      type_( (Type)params.get<int>( "type", (int)Type::Invalid ) ),
      FE( 0. ), FM( 0. ), GE( 0. ), GM( 0. )
    {
      if ( params.has<ParametersList>( "structureFunctions" ) )
        str_fun_ = strfun::StructureFunctionsHandler::get().build( params.get<ParametersList>( "structureFunctions" ) );
    }

    void
    Parameterisation::setStructureFunctions( const std::shared_ptr<strfun::Parameterisation>& sfmod )
    {
      str_fun_ = sfmod;
    }

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

    ArringtonEtAl::ArringtonEtAl( const ParametersList& params ) :
      mode_( params.get<int>( "mode" ) )
    {
      switch ( mode_ ) {
        case 0: // original
          a_e_ = { 3.439,-1.602, 0.068 };
          b_e_ = { 15.055,48.061,99.304,0.012,8.650 };
          a_m_ = { -1.465,1.260,0.262 };
          b_m_ = {  9.627,0.,0.,11.179,13.245 };
          break;
        case 1: // fit of quoted Ge+dGe values
          a_e_ = { 4.309,-1.108,-0.324 };
          b_e_ = { 15.340,58.321,124.11,3.927,0.589 };
          a_m_ = { -1.472,1.210,0.334 };
          b_m_ = {  9.486,0.,0., 9.440,15.416 };
          break;
        case 2: // fit of quted Ge-dGe values
          a_e_ = { 4.286,-1.281,-0.486 };
          b_e_ = { 16.308,54.535,138.03,7.005,0.014 };
          a_m_ = { -1.374,1.080,0.124 };
          b_m_ = { 10.003,0.,0., 7.680, 9.009 };
          break;
        case 3: // fit of quted Ge values
          a_e_ = { 4.109,-1.052,-0.375 };
          b_e_ = { 15.602,55.519,123.96,11.403,1.931 };
          a_m_ = { -1.436,1.196,0.210 };
          b_m_ = {  9.721,0.,0., 9.623,11.817 };
          break;
      }
    }

    void
    ArringtonEtAl::compute( double q2 )
    {
      const double tau = 0.25*q2/mp2_;

      double num_e = 1., den_e = 1.;
      for ( unsigned short i = 0; i < a_e_.size(); ++i )
        num_e += a_e_.at( i )*pow( tau, i+1 );
      for ( unsigned short i = 0; i < b_e_.size(); ++i )
        den_e += b_e_.at( i )*pow( tau, i+1 );
      GE = num_e/den_e;

      double num_m = 1., den_m = 1.;
      for ( unsigned short i = 0; i < a_m_.size(); ++i )
        num_m += a_m_.at( i )*pow( tau, i+1 );
      for ( unsigned short i = 0; i < b_m_.size(); ++i )
        den_m += b_m_.at( i )*pow( tau, i+1 );
      GM = num_m/den_m;
    }

    void
    BrashEtAl::compute( double q2 )
    {
      const double q = sqrt( q2 );
      GM = 1./( 1.+q*( 0.116+q*( 2.874+q*( 0.241+q*( 1.006+q*0.345 ) ) ) ) );

      const double r = std::min( 1., 1.-0.13*( q2-0.04 ) );
      if ( r < 0. ) {
        GM = GE = 0.;
        return;
      }
      GE = r*GM;
    }

    std::ostream&
    operator<<( std::ostream& os, const ff::Parameterisation& formfac )
    {
      return os << "FF[" << (int)formfac.type_ << "|" << (int)formfac.model_
        << "]{FE=" << formfac.FE << ",FM=" << formfac.FM << "}";
    }
  }
}
