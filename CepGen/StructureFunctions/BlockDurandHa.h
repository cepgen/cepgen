#ifndef CepGen_StructureFunctions_BlockDurandHa_h
#define CepGen_StructureFunctions_BlockDurandHa_h

#include "StructureFunctions.h"
#include <array>

namespace CepGen
{
  namespace SF
  {
    struct BlockDurandHaParameters
    {
      std::array<double,3> a, b;
      std::array<double,2> c;
      double n;
      /// Effective mass spread parameter
      double lambda;
      /// Asymptotic log-behaviour transition scale factor
      double mu2;
      /// Squared effective mass (~VM mass)
      double m2;
      static BlockDurandHaParameters standard();
    };
    StructureFunctions BlockDurandHa( double q2, double xbj, const BlockDurandHaParameters& params=BlockDurandHaParameters::standard() );
  }
}

#endif
