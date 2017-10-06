#include "Schaefer.h"

namespace CepGen
{
  namespace SF
  {
    Schaefer::Schaefer()
    {
#ifdef SchaeferF2
      luxlike_params_.amp = Constants::mp;
      luxlike_params_.am_pi = Constants::mpi;
      luxlike_params_.alpha_em = Constants::alphaEM;
      luxlike_params_.q2_cut = 9.;
      luxlike_params_.w2_hi = 4.;
      luxlike_params_.w2_lo = 3.;
      luxlike_params_.res_model = ChristyBosted;
      luxlike_params_.cont_model = GD11p;
#endif
    }

    StructureFunctions
    Schaefer::operator()( double q2, double xbj ) const
    {
      StructureFunctions luxlike;
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
