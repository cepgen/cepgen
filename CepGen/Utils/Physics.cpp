#include "CepGen/Utils/Physics.h"

namespace cepgen
{
  namespace utils
  {
    double mX2( double xbj, double q2, double mp2 )
    {
      return mp2+q2*( 1.-xbj )/xbj;
    }
  }
}
