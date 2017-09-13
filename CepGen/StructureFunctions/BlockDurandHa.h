#ifndef CepGen_StructureFunctions_BlockDurandHa_h
#define CepGen_StructureFunctions_BlockDurandHa_h

#include "StructureFunctions.h"

namespace CepGen
{
  namespace SF
  {
    struct BlockDurandHaParameters
    {
      std::vector<double> a, b, c;
      double xn, xlambda, mu2, xm2;
      static BlockDurandHaParameters standard();
    };
    StructureFunctions BlockDurandHa( double q2, double xbj, const BlockDurandHaParameters& params=BlockDurandHaParameters::standard() );
  }
}

#endif
