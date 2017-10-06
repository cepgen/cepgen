#ifndef CepGen_StructureFunctions_ALLM_h
#define CepGen_StructureFunctions_ALLM_h

#include "StructureFunctions.h"

namespace CepGen
{
  namespace SF
  {
    class ALLMParameterisation
    {
      private:
        struct Parameters {
          Parameters() :
            a( { 0., 0., 0. } ), b( { 0., 0., 0. } ), c( { 0., 0., 0. } ) {}
          Parameters( std::vector<double> c, std::vector<double> a, std::vector<double> b ) :
            a( a ), b( b ), c( c ) {}
          std::vector<double> a, b, c;
        };

      public:
        /// Pre-HERA data fit (694 data points)
        static ALLMParameterisation allm91();
        /// Fixed target and HERA photoproduction total cross sections (1356 points)
        static ALLMParameterisation allm97();
        static ALLMParameterisation hht_allm();
        static ALLMParameterisation hht_allm_ft();
        static ALLMParameterisation gd07p();
        static ALLMParameterisation gd11p();

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
    };

    StructureFunctions ALLM( double q2, double xbj, const ALLMParameterisation& param=ALLMParameterisation::allm97() );
  }
}

#endif
