#include "Schaefer.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/ParticleProperties.h"

namespace CepGen
{
  namespace SF
  {
    Schaefer::Schaefer()
    {
#ifdef SchaeferF2
      luxlike_params_.amp = ParticleProperties::mass( Proton );
      luxlike_params_.alpha_em = Constants::alphaEM;
      luxlike_params_.q2_cut = 9.;
      luxlike_params_.w2_hi = 4.;
      luxlike_params_.w2_lo = 3.;
      luxlike_params_.res_model = ChristyBosted;
      //luxlike_params_.cont_model = ALLM97;
      luxlike_params_.cont_model = GD11p;
      luxlike_params_.higher_twist = 0;
#endif
    }

    Schaefer
    Schaefer::operator()( double q2, double xbj ) const
    {
      Schaefer luxlike;
#ifndef SchaeferF2
      FatalError( "LUXlike structure functions cannot be computed "
                  "as the Fortran subroutine is not linked to this instance!" );
#else
      f2_fit_luxlike_( xbj, q2, luxlike.F2, luxlike.FL );
#endif
      return luxlike;
    }
  }
}
