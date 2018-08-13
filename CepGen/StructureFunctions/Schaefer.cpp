#include "CepGen/StructureFunctions/Schaefer.h"
#include "CepGen/StructureFunctions/StructureFunctionsBuilder.h"

#include "CepGen/Core/Exception.h"

#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/ParticleProperties.h"

namespace CepGen
{
  namespace SF
  {
    Schaefer::Parameterisation
    Schaefer::Parameterisation::standard()
    {
      Parameterisation par;
      par.q2_cut = 9.;
      par.w2_hi = 4.;
      par.w2_lo = 3.;
      par.resonances_model = StructureFunctionsBuilder::get( SF::Type::ChristyBosted );
      par.perturbative_model = StructureFunctionsBuilder::get( SF::Type::MSTWgrid );
      par.continuum_model = StructureFunctionsBuilder::get( SF::Type::GD11p );
      par.higher_twist = 0;
      return par;
    }

    Schaefer::Parameterisation
    Schaefer::Parameterisation::cteq()
    {
      Parameterisation par;
      par.q2_cut = 9.;
      par.w2_hi = 4.;
      par.w2_lo = 3.;
      par.resonances_model = StructureFunctionsBuilder::get( SF::Type::ChristyBosted );
      par.perturbative_model = StructureFunctionsBuilder::get( SF::Type::MSTWgrid );
      par.continuum_model = StructureFunctionsBuilder::get( SF::Type::GD11p );
      par.higher_twist = 0;
      return par;
    }

    Schaefer::Schaefer( const Schaefer::Parameterisation& params ) :
      StructureFunctions( Type::Schaefer ),
      params( params ), initialised_( false ), inv_omega_range_( -1. )
    {}

    void
    Schaefer::initialise()
    {
      CG_INFO( "LUXlike" ) << "LUXlike structure functions evaluator successfully initialised.\n"
        << " * Q² cut:             " << params.q2_cut << " GeV²\n"
        << " * W² ranges:          " << params.w2_lo << " GeV² / " << params.w2_hi << " GeV²\n"
        << " * resonance model:    " << params.resonances_model->type << "\n"
        << " * perturbative model: " << params.perturbative_model->type << "\n"
        << " * continuum model:    " << params.continuum_model->type << "\n"
        << " * higher-twist?       " << std::boolalpha << params.higher_twist;
      inv_omega_range_ = 1./( params.w2_hi-params.w2_lo );
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

      StructureFunctions sel_sf;
      if ( q2 < params.q2_cut ) {
        if ( w2 < params.w2_lo )
          sel_sf = ( *params.resonances_model )( xbj, q2 );
        else if ( w2 < params.w2_hi ) {
          auto sf_r = ( *params.resonances_model )( xbj, q2 );
          auto sf_c = ( *params.continuum_model )( xbj, q2 );
          sf_r.computeFL( xbj, q2 );
          sf_c.computeFL( xbj, q2 );
          const double r = rho( w2 );
          F2 = r*sf_c.F2 + ( 1.-r )*sf_r.F2;
          FL = r*sf_c.FL + ( 1.-r )*sf_r.FL;
          return *this;
        }
        else
          sel_sf = ( *params.continuum_model )( xbj, q2 );
      }
      else {
        if ( w2 < params.w2_hi )
          sel_sf = ( *params.continuum_model )( xbj, q2 );
        else {
          auto sf_p = ( *params.perturbative_model )( xbj, q2 );
          F2 = sf_p.F2;
          sf_p.computeFL( xbj, q2 );
          FL = sel_sf.FL;
          if ( params.higher_twist )
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
          << params.w2_lo << " / " << params.w2_hi << " GeV²!";
      const double omega = ( w2-params.w2_lo )*inv_omega_range_;
      const double omega2 = omega*omega;
      return 2.*omega2-omega*omega;
    }
  }
}
