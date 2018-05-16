#include "CepGen/Core/GridParameters.h"
#include "CepGen/Core/Exception.h"

namespace CepGen
{
  const unsigned short GridParameters::max_dimensions_ = 15;
  const unsigned short GridParameters::mbin_ = 3;
  const double GridParameters::inv_mbin_ = 1./mbin_;

  GridParameters::GridParameters() :
    gen_prepared( false ), f_max_global( 0. ), f_max_diff( 0. )
  {}
}
