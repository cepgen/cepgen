#ifndef CepGen_StructureFunctions_FioreBrasse_h
#define CepGen_StructureFunctions_FioreBrasse_h

#include "CepGen/StructureFunctions/StructureFunctions.h"

#include <vector>

namespace CepGen
{
  namespace SF
  {
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
        /// Fiore-Brasse proton structure functions (F.W Brasse et al., DESY 76/11 (1976),
        /// http://dx.doi.org/10.1016/0550-3213(76)90231-5)
        explicit FioreBrasse( const FioreBrasse::Parameterisation& params = FioreBrasse::Parameterisation::standard() );
        /// \param[in] q2 Squared 4-momentum transfer
        /// \param[in] xbj Bjorken's x
        /// \cite Brasse1976413
        FioreBrasse& operator()( double xbj, double q2 ) override;
        FioreBrasse& operator()( double xbj, double q2, bool old );

        double W1, W2;
        Parameterisation params;
    };
  }
}

#endif
