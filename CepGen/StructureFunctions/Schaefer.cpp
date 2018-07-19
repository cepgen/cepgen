#include "CepGen/StructureFunctions/Schaefer.h"
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
      par.amp = mp_;
      par.alpha_em = Constants::alphaEM;
      par.q2_cut = 9.;
      par.w2_hi = 4.;
      par.w2_lo = 3.;
      par.res_model = (int)SF::Type::ChristyBosted;
      par.cont_model = (int)SF::Type::GD11p;
      par.higher_twist = 0;
      return par;
    }

    Schaefer::Schaefer( const Schaefer::Parameterisation& params ) :
      StructureFunctions( Type::Schaefer ),
      params( params ), initialised_( false )
    {}

    void
    Schaefer::initialise()
    {
#ifdef SchaeferF2
      luxlike_params_ = params;
      CG_INFO( "Schaefer" ) << "LUXlike structure functions evaluator successfully initialised.\n"
        << " * proton mass:     " << params.amp << " GeV/c²\n"
        << " * alpha(em):       " << params.alpha_em << "\n"
        << " * Q² cut:          " << params.q2_cut << " GeV²\n"
        << " * W² ranges:       " << params.w2_lo << " GeV² / " << params.w2_hi << " GeV²\n"
        << " * resonance model: " << (SF::Type)params.res_model << "\n"
        << " * continuum model: " << (SF::Type)params.cont_model << "\n"
        << " * higher-twist?    " << std::boolalpha << (bool)params.higher_twist;
      initialised_ = true;
#endif
    }

    Schaefer&
    Schaefer::operator()( double q2, double xbj )
    {
#ifndef SchaeferF2
      throw CG_FATAL( "Schaefer" )
        << "LUXlike structure functions cannot be computed "
        << "as the Fortran subroutine is not linked to this instance!";
#else
      if ( !initialised_ )
        initialise();
      std::pair<double,double> nv = { q2, xbj };
      if ( nv == old_vals_ )
        return *this;
      old_vals_ = nv;

      f2_fit_luxlike_( xbj, q2, F2, FL );
#endif
      return *this;
    }
  }
}
