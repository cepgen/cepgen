#ifndef CepGen_StructureFunctions_BlockDurandHa_h
#define CepGen_StructureFunctions_BlockDurandHa_h

#include "StructureFunctions.h"
#include <array>

namespace CepGen
{
  namespace SF
  {
    class BlockDurandHa : public StructureFunctions
    {
      public:
        struct Parameterisation
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
          static Parameterisation standard();
        };
        explicit BlockDurandHa( const BlockDurandHa::Parameterisation& params = BlockDurandHa::Parameterisation::standard() );
        BlockDurandHa& operator()( double q2, double xbj ) override;

      private:
        Parameterisation params_;
    };
  }
}

#endif
