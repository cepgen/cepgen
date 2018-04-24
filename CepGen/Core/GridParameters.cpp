#include "CepGen/Core/GridParameters.h"
#include "CepGen/Core/Exception.h"

namespace CepGen
{
  GridParameters::GridParameters() :
    gen_prepared( false ), f_max_global( 0. ), f_max_diff( 0. )
  {}
}
