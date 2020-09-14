#include "CepGen/Physics/FormFactors.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/Constants.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/String.h"

#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/StructureFunctions/SuriYennie.h"

#include <cmath>
#include <cassert>
//#include <iosfwd>

namespace cepgen
{
  namespace ff
  {
    Parameterisation::Parameterisation() :
      NamedModule<int>( ParametersList() ),
      type_( Type::Invalid ),
      mp_( PDG::get().mass( PDG::proton ) ), mp2_( mp_*mp_ ),
      last_q2_( -1. ),
      FE( 0. ), FM( 0. ), GE( 0. ), GM( 0. )
    {}

    Parameterisation::Parameterisation( const Parameterisation& param ) :
      NamedModule<int>( param.parameters() ),
      type_( param.type_ ),
      mp_( param.mp_ ), mp2_( param.mp2_ ),
      last_q2_( -1. ),
      FE( param.FE ), FM( param.FM ), GE( param.GE ), GM( param.GM )
    {}

    Parameterisation::Parameterisation( const ParametersList& params ) :
      NamedModule<int>( params ),
      type_( (Type)params.get<int>( "type", (int)Type::Invalid ) ),
      mp_( PDG::get().mass( PDG::proton ) ), mp2_( mp_*mp_ ),
      last_q2_( -1. ),
      FE( 0. ), FM( 0. ), GE( 0. ), GM( 0. )
    {}

    void
    Parameterisation::setStructureFunctions( strfun::Parameterisation* sfmod )
    {
      str_fun_.reset( sfmod );
    }

    double
    Parameterisation::tau( double q2 ) const
    {
      if ( mp2_ <= 0. )
        throw CG_FATAL( "FormFactors:tau" )
          << "Invalid proton mass! check the form factors constructor!";
      return 0.25*q2/mp2_;
    }

    Parameterisation&
    Parameterisation::operator()( double q2, double mi2, double mf2 )
    {
      last_q2_ = q2;
      switch ( type_ ) {
        case Type::Invalid:
        case Type::CompositeScalar:
          throw CG_FATAL( "FormFactors" )
            << type_ << " mode is not yet supported!";
        case Type::PointLikeScalar:
          FE = 1., FM = 0.; break;
        case Type::PointLikeFermion:
          FE = FM = 1.; break;
        case Type::ProtonElastic: {
          compute( q2 );
          const double GE2 = GE*GE, GM2 = GM*GM;
          FE = ( 4.*mp2_*GE2+q2*GM2 ) / ( 4.*mp2_+q2 );
          FM = GM2;
        } break;
        case Type::ProtonInelastic: {
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

    class StandardDipole : public Parameterisation
    {
      public:
        using Parameterisation::Parameterisation;
        static std::string description() { return "Standard dipole"; }

      private:
        void compute( double q2 ) override {
          GE = pow( 1.+q2/0.71, -2. );
          GM = MU*GE;
        }
    };

    //------------------------------------------------------------------

    class ArringtonEtAl : public Parameterisation
    {
      public:
        ArringtonEtAl( const ParametersList& );
        static std::string description() { return "Arrington et al."; }

      private:
        void compute( double q2 ) override;
        const int mode_;
        std::vector<double> a_e_, b_e_;
        std::vector<double> a_m_, b_m_;
    };

    ArringtonEtAl::ArringtonEtAl( const ParametersList& params ) :
      Parameterisation( params ),
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
      const double tau_val = tau( q2 );

      double num_e = 1., den_e = 1.;
      for ( unsigned short i = 0; i < a_e_.size(); ++i )
        num_e += a_e_.at( i )*pow( tau_val, i+1 );
      for ( unsigned short i = 0; i < b_e_.size(); ++i )
        den_e += b_e_.at( i )*pow( tau_val, i+1 );
      GE = num_e/den_e;

      double num_m = 1., den_m = 1.;
      for ( unsigned short i = 0; i < a_m_.size(); ++i )
        num_m += a_m_.at( i )*pow( tau_val, i+1 );
      for ( unsigned short i = 0; i < b_m_.size(); ++i )
        den_m += b_m_.at( i )*pow( tau_val, i+1 );
      GM = MU*num_m/den_m;
    }

    //------------------------------------------------------------------

    class BrashEtAl : public Parameterisation
    {
      public:
        using Parameterisation::Parameterisation;
        static std::string description() { return "Brash et al."; }

      private:
        static constexpr float MAX_Q2 = 7.7;
        void compute( double q2 ) override {
          if ( q2 > MAX_Q2 )
            CG_WARNING( "BrashEtAl" )
            << "Q² = " << q2 << " > " << MAX_Q2 << " GeV² = max(Q²).\n\t"
            << "Brash et al. FF parameterisation not designed for high-Q² values.";
          const double q = sqrt( q2 );
          GM = 1./( 1.+q*( 0.116+q*( 2.874+q*( 0.241+q*( 1.006+q*0.345 ) ) ) ) );

          const double r = std::min( 1., 1.-0.13*( q2-0.04 ) );
          if ( r < 0. ) {
            GM = GE = 0.;
            return;
          }
          GE = r*GM;
          GM *= MU;
        }
    };

    //------------------------------------------------------------------

    class MergellEtAl : public Parameterisation
    {
      public:
        MergellEtAl( const ParametersList& );
        static std::string description() { return "Mergell et al."; }

      private:
        void compute( double q2 ) override;
        static constexpr double Q2_RESCL = 9.733, INV_DENUM = 1./0.350;
        static constexpr double EXPO = 2.148;
        const std::vector<double> par1_, par2_;
    };

    MergellEtAl::MergellEtAl( const ParametersList& params ) :
      Parameterisation( params ),
      par1_( params.get<std::vector<double> >( "par1", { 1.0317, 0.0875, 0.3176, 0.5496 } ) ),
      par2_( params.get<std::vector<double> >( "par2", { 5.7824, 0.3907, 0.1422, 0.5362 } ) )
    {
      assert( par1_.size() == 4 );
      assert( par2_.size() == 4 );
    }

    void
    MergellEtAl::compute( double q2 )
    {
      const double log1 = std::pow( log( ( Q2_RESCL+q2 )*INV_DENUM ), -EXPO );
      const double d1_1 = 0.611+q2, d2_1 = 1.039+q2, d3_1 = 2.560+q2;

      const double Fs1 = (  9.464/d1_1-9.054/d2_1-0.410/d3_1 )*log1;
      const double Fs2 = ( -1.549/d1_1+1.985/d2_1-0.436/d3_1 )*log1;

      const double log2 = std::pow( log( ( Q2_RESCL-0.500 )*INV_DENUM ), +EXPO );
      const double log3 = std::pow( log( ( Q2_RESCL-0.400 )*INV_DENUM ), +EXPO );
      const double d1_2 = 2.103+q2, d2_2 = 2.734+q2, d3_2 = 2.835+q2;

      const double Fv1= ( 0.5*( par1_.at(0)*log2+par1_.at(1)*log3*std::pow( 1.+q2/par1_.at(2), -2 ) )/( 1.+q2/par1_.at(3) )-38.885/d1_2+425.007/d2_2-389.742/d3_2 )*log1;
      const double Fv2= ( 0.5*( par2_.at(0)*log2+par2_.at(1)*log3/        ( 1.+q2/par2_.at(2)     ) )/( 1.+q2/par2_.at(3) )-73.535/d1_2+ 83.211/d2_2- 29.467/d3_2 )*log1;

      const double F1 = Fv1+Fs1, F2 = Fv2+Fs2;

      GE = F1-tau( q2 )*F2;
      GM = F1+F2;
    }

    //------------------------------------------------------------------

    std::ostream&
    operator<<( std::ostream& os, const Parameterisation& formfac )
    {
      os << formfac.description();
      if ( formfac.last_q2_ >= 0. )
        os << "(Q²=" << formfac.last_q2_ << " GeV²): "
           << "FE=" << formfac.FE << ",FM=" << formfac.FM;
      return os;
    }

    std::ostream&
    operator<<( std::ostream& os, const Type& type )
    {
      switch ( type ) {
        case Type::Invalid: return os << "{invalid}";
        case Type::ProtonElastic: return os << "el.proton";
        case Type::PointLikeScalar: return os << "gen.scalar";
        case Type::PointLikeFermion: return os << "gen.fermion";
        case Type::CompositeScalar: return os << "comp.scalar";
        case Type::ProtonInelastic: return os << "inel.proton";
      }
      return os;
    }

    std::ostream&
    operator<<( std::ostream& os, const Model& mod )
    {
      switch ( mod ) {
        case Model::Invalid: return os << "{invalid}";
        case Model::StandardDipole: return os << "std.dipole";
        case Model::ArringtonEtAl: return os << "Arrington etc.";
        case Model::BrashEtAl: return os << "Brash etc.";
        case Model::MergellEtAl: return os << "Mergell etc.";
      }
      return os;
    }
  }
}

REGISTER_FF_MODEL( StandardDipole, ff::StandardDipole )
REGISTER_FF_MODEL( ArringtonEtAl, ff::ArringtonEtAl )
REGISTER_FF_MODEL( BrashEtAl, ff::BrashEtAl )
REGISTER_FF_MODEL( MergellEtAl, ff::MergellEtAl )
