#ifndef CepGen_StructureFunctions_FioreBrasse_h
#define CepGen_StructureFunctions_FioreBrasse_h

#include "CepGen/StructureFunctions/StructureFunctions.h"

#include <vector>

namespace CepGen
{
  namespace SF
  {
    ///\f${\cal W}_{1,2}\f$ structure functions parameterisation by Fiore et al \cite Fiore:2002re and Brasse et al \cite Brasse:1976bf
    class FioreBrasse : public StructureFunctions
    {
      public:
        struct Parameterisation
        {
          static Parameterisation standard();
          static Parameterisation alternative();

          struct ResonanceParameters {
//            ResonanceParameters( double a0, double a1, double a2, double a, double q02, float spin ) :
//              alpha0( a0 ), alpha1( a1 ), alpha2( a2 ), a( a ), q02( q02 ), spin( spin ) {}
            double alpha0, alpha1, alpha2, a, q02;
            float spin;
          };

          std::vector<ResonanceParameters> resonances;
          double s0, norm;
        };
        /// Fiore \cite Fiore:2002re and Brasse \cite Brasse:1976bf proton structure functions
        explicit FioreBrasse( const FioreBrasse::Parameterisation& params = FioreBrasse::Parameterisation::standard() );
        FioreBrasse& operator()( double xbj, double q2 ) override;
        /// Old implementation from LPAIR
        FioreBrasse& operator()( double xbj, double q2, bool old );

        double W1, W2;
        Parameterisation params;
    };
  }
}

#endif
