#ifndef CepGen_StructureFunctions_BlockDurandHa_h
#define CepGen_StructureFunctions_BlockDurandHa_h

#include "CepGen/StructureFunctions/StructureFunctions.h"
#include <vector>

namespace cepgen
{
  namespace strfun
  {
    /// \f$F_2\f$ parameterisation from Block, Durand, and Ha \cite Block:2014kza
    class BlockDurandHa : public Parameterisation
    {
      public:
        explicit BlockDurandHa( const ParametersList& params = ParametersList() );
        BlockDurandHa& operator()( double xbj, double q2 ) override;

      private:
        std::vector<double> a_, b_, c_;
        double n_;
        /// Effective mass spread parameter
        double lambda_;
        /// Asymptotic log-behaviour transition scale factor
        double mu2_;
        /// Squared effective mass (~VM mass)
        double m2_;
    };
  }
}

#endif
