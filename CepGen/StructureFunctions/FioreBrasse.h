#ifndef CepGen_StructureFunctions_FioreBrasse_h
#define CepGen_StructureFunctions_FioreBrasse_h

#include "CepGen/StructureFunctions/StructureFunctions.h"

#include <vector>

namespace cepgen
{
  namespace strfun
  {
    ///\f${\cal W}_{1,2}\f$ structure functions parameterisation by Fiore et al \cite Fiore:2002re and Brasse et al \cite Brasse:1976bf
    class FioreBrasse : public Parameterisation
    {
      public:
        /// General parameters for this modelling
        struct Parameters
        {
          static Parameters standard();
          static Parameters alternative();
          /// Description of a single resonance in the modelling
          struct Resonance {
            double alpha0, alpha1, alpha2, a, q02;
            float spin;
          };
          /// All resonances considered in this modelling
          std::vector<Resonance> resonances;
          double s0, norm;
        };
        /// Fiore \cite Fiore:2002re and Brasse \cite Brasse:1976bf proton structure functions
        explicit FioreBrasse( const ParametersList& params = ParametersList() );
        FioreBrasse& operator()( double xbj, double q2 ) override;

      private:
        Parameters params_;
    };
  }
}

#endif
