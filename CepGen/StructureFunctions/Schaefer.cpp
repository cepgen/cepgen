#include "Schaefer.h"
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
      par.res_model = ChristyBosted;
      par.cont_model = GD11p;
      par.higher_twist = 0;
      return par;
    }

    Schaefer::Schaefer( const Schaefer::Parameterisation& params )
    {
#ifdef SchaeferF2
      luxlike_params_ = params;
#endif
    }

    Schaefer
    Schaefer::operator()( double q2, double xbj ) const
    {
      Schaefer luxlike;
#ifndef SchaeferF2
      throw FatalError( "Schaefer" )
        << "LUXlike structure functions cannot be computed "
        << "as the Fortran subroutine is not linked to this instance!";
#else
      f2_fit_luxlike_( xbj, q2, luxlike.F2, luxlike.FL );
#endif
      return luxlike;
    }
  }
}
