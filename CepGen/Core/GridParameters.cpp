#include "CepGen/Core/GridParameters.h"
#include "CepGen/Core/Exception.h"
#include <cmath> // pow

namespace cepgen
{
  GridParameters::GridParameters( unsigned short ndim ) :
    max( pow( M_BIN, ndim ) ), gen_prepared( false ),
    f_max( max, 0. ), f_max_global( 0. ), f_max_diff( 0. ),
    num( max ), r_boxes( 0 )
  {
    if ( ndim > MAX_DIM )
      throw CG_FATAL( "GridParameters" ) << "Phase space too large!\n\t"
        << "Either reduce the number of integration dimensions, or\n\t"
        << "increase the GridParameters::MAX_DIM parameter (not recommended).";
  }
}
