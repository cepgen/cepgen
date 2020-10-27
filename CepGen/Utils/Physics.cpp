#include "CepGen/Utils/Physics.h"

namespace cepgen
{
  namespace utils
  {
    double mX2( double xbj, double q2, double mp2 )
    {
      return mp2+q2*( 1.-xbj )/xbj;
    }

    double xBj( double q2, double mp2, double mx2 )
    {
      return q2/( q2-mp2+mx2 );
    }
  }
}
