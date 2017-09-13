#ifndef CepGen_Physics_StructureFunctions_h
#define CepGen_Physics_StructureFunctions_h

#include "Particle.h"
#include <complex>

extern "C"
{
  extern void grv95lo_( float&, float&, float&, float&, float&, float&, float&, float& );
  //extern void grv95lo_( double&, double&, double&, double&, double&, double&, double&, double& );
}

namespace CepGen
{
  /// Proton structure function to be used in the outgoing state description
  /// \note Values correspond to the LPAIR legacy steering card values
  enum StructureFunctionsType {
    Electron = 1,
    ElasticProton = 2,
    SuriYennie = 11,
    SuriYennieLowQ2 = 12,
    SzczurekUleshchenko = 15,
    FioreVal = 101,
    FioreSea = 102,
    Fiore = 103,
    ALLM
  };
  /// Human-readable format of a structure function type
  std::ostream& operator<<( std::ostream& os, const StructureFunctionsType& sf );

  class StructureFunctions
  {
    public:
      StructureFunctions( double f1=0.0, double f2=0.0 ) : F1( f1 ), F2( f2 ) {}

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
      static StructureFunctions FioreBrasse( double q2, double xbj );
      static StructureFunctions FioreBrasseOld( double q2, double xbj );
      struct SuriYennieParameterisation {
        double C1, C2, D1, rho2, Cp, Bp;
        // values extracted from experimental fits
        static SuriYennieParameterisation standard() {
          SuriYennieParameterisation p; p.C1 = 0.86926; p.C2 = 2.23422; p.D1 = 0.12549; p.rho2 = 0.585; p.Cp = 0.96; p.Bp = 0.63; return p;
        }
        static SuriYennieParameterisation alternative() {
          SuriYennieParameterisation p; p.C1 = 0.6303; p.C2 = 2.3049; p.D1 = 0.04681; p.rho2 = 1.05; p.Cp = 1.23; p.Bp = 0.61; return p;
        }
      };
      static StructureFunctions SuriYennie( double q2, double xbj, const SuriYennieParameterisation& param=SuriYennieParameterisation::standard() );
      static StructureFunctions SzczurekUleshchenko( double q2, double xbj );

      struct ALLMParameterisation {
        struct Parameters {
          Parameters() :
            a1( 0. ), a2( 0. ), a3( 0. ), b1( 0. ), b2( 0. ), b3( 0. ), c1( 0. ), c2( 0. ), c3( 0. ) {}
          Parameters( double c1, double c2, double c3, double a1, double a2, double a3, double b1, double b2, double b3 ) :
            a1( a1 ), a2( a2 ), a3( a3 ), b1( b1 ), b2( b2 ), b3( b3 ), c1( c1 ), c2( c2 ), c3( c3 ) {}
          double a1, a2, a3, b1, b2, b3, c1, c2, c3;
        };
        Parameters pomeron, reggeon;
        double m02, mp2, mr2;
        double q02, lam2;
        static ALLMParameterisation standard() {
          ALLMParameterisation p;
          p.pomeron = Parameters( 0.28067, 0.22291, 2.1979, -0.0808, -0.44812, 1.1709, 0.36292, 1.8917, 1.8439 );
          p.reggeon = Parameters( 0.80107, 0.97307, 3.4924, 0.58400, 0.37888, 2.6063, 0.01147, 3.7582, 0.49338 );
          p.m02 = 0.31985;
          p.mp2 = 49.457;
          p.mr2 = 0.15052;
          p.q02 = 0.52544;
          p.lam2 = 0.06526;
          return p;
        }
      };
      static StructureFunctions ALLM( double q2, double xbj, const ALLMParameterisation& param=ALLMParameterisation::standard() );
      double F1, F2;
      double FM;
  };
  /// Human-readable format of a structure function object
  std::ostream& operator<<( std::ostream& os, const StructureFunctions& sf );
}

#endif
