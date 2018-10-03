#ifndef CepGen_StructureFunctions_ALLM_h
#define CepGen_StructureFunctions_ALLM_h

#include "StructureFunctions.h"
#include "SigmaRatio.h"

#include <vector>

namespace CepGen
{
  namespace SF
  {
    /// \f$F_{2/L}\f$ parameterisation by Abramowicz, Levin, Levy, and Maor \cite Abramowicz:1991xz\cite Abramowicz:1997ms
    class ALLM : public StructureFunctions
    {
      public:
        class Parameterisation
        {
          private:
            struct Parameters {
              Parameters() :
                a( { 0., 0., 0. } ), b( { 0., 0., 0. } ), c( { 0., 0., 0. } ) {}
              Parameters( const std::vector<double>& c, const std::vector<double>& a, const std::vector<double>& b ) :
                a( a ), b( b ), c( c ) {}
              std::vector<double> a, b, c;
            };

          public:
            Parameterisation() :
              m02( 0. ), mp2( 0. ), mr2( 0. ), q02( 0. ), lambda2( 0. ) {}
            /// Pre-HERA data fit (694 data points)
            static Parameterisation allm91();
            /// Fixed target and HERA photoproduction total cross sections (1356 points)
            static Parameterisation allm97();
            static Parameterisation hht_allm();
            static Parameterisation hht_allm_ft();
            static Parameterisation gd07p();
            static Parameterisation gd11p();

            Parameters pomeron, reggeon;
            /// Effective photon squared mass
            double m02;
            /// Effective pomeron squared mass
            double mp2;
            /// Effective reggeon squared mass
            double mr2;
            double q02;
            /// Squared QCD scale
            double lambda2;
            Type type;
        };

        explicit ALLM( const ALLM::Parameterisation& param = ALLM::Parameterisation::allm97() );
        ALLM& operator()( double xbj, double q2 ) override;

      private:
        Parameterisation params_;
    };
  }
}

#endif
