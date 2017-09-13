#ifndef CepGen_StructureFunctions_FioreBrasse_h
#define CepGen_StructureFunctions_FioreBrasse_h

#include "StructureFunctions.h"
#include <complex>

namespace CepGen
{
  namespace SF
  {
    struct FioreBrasseParameterisation {
      struct ResonanceParameters {
        ResonanceParameters( double a0, double a1, double a2, double a, double q02, float spin, bool on=true ) :
          alpha0( a0 ), alpha1( a1 ), alpha2( a2 ), a( a ), q02( q02 ), spin( spin ), enabled( on ) {}
        double alpha0, alpha1, alpha2, a, q02;
        float spin;
        bool enabled;
      };
      std::vector<ResonanceParameters> resonances;
      double s0, norm;

      static FioreBrasseParameterisation standard() {
        FioreBrasseParameterisation p;
        p.s0 = 1.14;
        p.norm = 0.021;
        p.resonances.emplace_back( -0.8377, 0.95, 0.1473, 1.0, 2.4617, 3./2. );
        p.resonances.emplace_back( -0.37, 0.95, 0.1471, 0.5399, 2.4617, 5./2. );
        p.resonances.emplace_back( 0.0038, 0.85, 0.1969, 4.2225, 1.5722, 3./2. );
        p.resonances.emplace_back( 0.5645, 0.1126, 1.3086, 19.2694, 4.5259, 1. );
        return p;
      }
      static FioreBrasseParameterisation alternative() {
        FioreBrasseParameterisation p;
        p.s0 = 1.2871;
        p.norm = 0.0207;
        p.resonances.emplace_back( -0.8070, 0.9632, 0.1387, 1.0, 2.6066, 3./2. );
        p.resonances.emplace_back( -0.3640, 0.9531, 0.1239, 0.6086, 2.6066, 5./2. );
        p.resonances.emplace_back( -0.0065, 0.8355, 0.2320, 4.7279, 1.4828, 3./2. );
        p.resonances.emplace_back( 0.5484, 0.1373, 1.3139, 14.7267, 4.6041, 1. );
        return p;
      }
    };
    /// Fiore-Brasse proton structure functions (F.W Brasse et al., DESY 76/11 (1976),
    /// http://dx.doi.org/10.1016/0550-3213(76)90231-5)
    /// \param[in] q2 Squared 4-momentum transfer
    /// \param[in] xbj Bjorken's x
    /// \cite Brasse1976413
    StructureFunctions FioreBrasse( double q2, double xbj );
    StructureFunctions FioreBrasseOld( double q2, double xbj );
  }
}

#endif
