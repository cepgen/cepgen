#ifndef CepGen_StructureFunctions_FioreBrasse_h
#define CepGen_StructureFunctions_FioreBrasse_h

#include "CepGen/Physics/Constants.h"
#include "CepGen/StructureFunctions/StructureFunctions.h"

#include <complex>
#include <vector>

namespace CepGen
{
  namespace SF
  {
    class FioreBrasse
    {
      public:
        struct Parameterisation
        {
          static Parameterisation standard();
          static Parameterisation alternative();

          struct ResonanceParameters {
            ResonanceParameters( double a0, double a1, double a2, double a, double q02, float spin, bool on=true ) :
              alpha0( a0 ), alpha1( a1 ), alpha2( a2 ), a( a ), q02( q02 ), spin( spin ), enabled( on ) {}
            double alpha0, alpha1, alpha2, a, q02;
            float spin;
            bool enabled;
          };

          std::vector<ResonanceParameters> resonances;
          double s0, norm;
        };
        /// Fiore-Brasse proton structure functions (F.W Brasse et al., DESY 76/11 (1976),
        /// http://dx.doi.org/10.1016/0550-3213(76)90231-5)
        FioreBrasse( const FioreBrasse::Parameterisation& params = FioreBrasse::Parameterisation::standard() ) : params_( params ) {}
        /// \param[in] q2 Squared 4-momentum transfer
        /// \param[in] xbj Bjorken's x
        /// \cite Brasse1976413
        StructureFunctions operator()( double q2, double xbj ) const;
        StructureFunctions operator()( double q2, double xbj, bool old ) const;

      private:
        Parameterisation params_;
    };
  }
}

#endif
