#ifndef CepGen_StructureFunctions_BlockDurandHa_h
#define CepGen_StructureFunctions_BlockDurandHa_h

#include "CepGen/StructureFunctions/StructureFunctions.h"
#include <array>

namespace CepGen
{
  namespace SF
  {
    /// \f$F_2\f$ parameterisation from Block, Durand, and Ha \cite Block:2014kza
    class BlockDurandHa : public Parameterisation
    {
      public:
        struct Parameters
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
          static Parameters standard();
        };
        explicit BlockDurandHa( const Parameters& params = Parameters::standard() );
        BlockDurandHa& operator()( double xbj, double q2 ) override;

      private:
        Parameters params_;
    };
  }
}

#endif
