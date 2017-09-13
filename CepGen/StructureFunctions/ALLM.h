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
        static ALLMParameterisation allm91();
        static ALLMParameterisation allm97();

        Parameters pomeron, reggeon;
        double m02, mp2, mr2;
        double q02, lam2;
    };

    StructureFunctions ALLM( double q2, double xbj, const ALLMParameterisation& param=ALLMParameterisation::allm97() );
  }
}

#endif
