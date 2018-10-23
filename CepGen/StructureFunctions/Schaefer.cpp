#include "CepGen/StructureFunctions/Schaefer.h"
#include "CepGen/StructureFunctions/Partonic.h"

#include "CepGen/Core/Exception.h"

#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/ParticleProperties.h"

namespace cepgen
{
  namespace strfun
  {
    Schaefer::Schaefer( const ParametersList& params ) :
      Parameterisation( params ),
      q2_cut_( params.get<double>( "Q2cut", 9. ) ),
      w2_lim_( params.get<std::vector<double> >( "W2limits", { 3., 4. } ) ),
      higher_twist_( params.get<bool>( "higherTwist", true ) ),
      resonances_model_  ( Parameterisation::build( params.get<ParametersList>( "resonancesSF", ParametersList().set<int>( "id", (int)Type::ChristyBosted ) ) ) ),
      perturbative_model_( Parameterisation::build( params.get<ParametersList>( "perturbativeSF", ParametersList().set<int>( "id", (int)Type::MSTWgrid ) ) ) ),
      continuum_model_   ( Parameterisation::build( params.get<ParametersList>( "continuumSF", ParametersList().set<int>( "id", (int)Type::GD11p ) ) ) ),
      initialised_( false ), inv_omega_range_( -1. )
    {}

    std::string
    Schaefer::description() const
    {
      std::ostringstream os;
      os << "LUXlike{"
         << "r=" << *resonances_model_ << ","
         << "p=" << *perturbative_model_ << ","
         << "c=" << *continuum_model_;
      if ( higher_twist_ )
        os << ",HT";
      os << "}";
      return os.str();
    }

    void
    Schaefer::initialise()
    {
      CG_DEBUG( "LUXlike" ) << "LUXlike structure functions evaluator successfully initialised.\n"
        << " * Q² cut:             " << q2_cut_ << " GeV²\n"
        << " * W² ranges:          " << w2_lim_.at( 0 ) << " GeV² / " << w2_lim_.at( 1 ) << " GeV²\n"
        << " * resonance model:    " << *resonances_model_ << "\n"
        << " * perturbative model: " << *perturbative_model_ << "\n"
        << " * continuum model:    " << *continuum_model_ << "\n"
        << " * higher-twist?       " << std::boolalpha << higher_twist_;
      inv_omega_range_ = 1./( w2_lim_.at( 1 )-w2_lim_.at( 0 ) );
      initialised_ = true;
    }

    Schaefer&
    Schaefer::operator()( double xbj, double q2 )
    {
      if ( !initialised_ )
        initialise();

      std::pair<double,double> nv = { xbj, q2 };
      if ( nv == old_vals_ )
        return *this;
      old_vals_ = nv;

      const double w2 = mp2_+q2*( 1.-xbj )/xbj;

      strfun::Parameterisation sel_sf;
      if ( q2 < q2_cut_ ) {
        if ( w2 < w2_lim_.at( 0 ) )
          sel_sf = ( *resonances_model_ )( xbj, q2 );
        else if ( w2 < w2_lim_.at( 1 ) ) {
          auto sf_r = ( *resonances_model_ )( xbj, q2 );
          auto sf_c = ( *continuum_model_ )( xbj, q2 );
          sf_r.computeFL( xbj, q2 );
          sf_c.computeFL( xbj, q2 );
          const double r = rho( w2 );
          F2 = r*sf_c.F2 + ( 1.-r )*sf_r.F2;
          FL = r*sf_c.FL + ( 1.-r )*sf_r.FL;
          return *this;
        }
        else
          sel_sf = ( *continuum_model_ )( xbj, q2 );
      }
      else {
        if ( w2 < w2_lim_.at( 1 ) )
          sel_sf = ( *continuum_model_ )( xbj, q2 );
        else {
          auto sf_p = ( *perturbative_model_ )( xbj, q2 );
          F2 = sf_p.F2;
          sf_p.computeFL( xbj, q2 );
          FL = sf_p.FL;
          if ( higher_twist_ )
            F2 *= ( 1.+5.5/q2 );
          return *this;
        }
      }

      F2 = sel_sf( xbj, q2 ).F2;
      sel_sf.computeFL( xbj, q2 );
      FL = sel_sf.FL;

      return *this;
    }

    double
    Schaefer::rho( double w2 ) const
    {
      if ( inv_omega_range_ <= 0. )
        throw CG_FATAL( "LUXlike" ) << "Invalid W² limits: "
          << w2_lim_.at( 0 ) << " / " << w2_lim_.at( 1 ) << " GeV²!";
      const double omega = ( w2-w2_lim_.at( 0 ) )*inv_omega_range_;
      const double omega2 = omega*omega;
      return 2.*omega2-omega*omega;
    }
  }
}
