#include "CepGen/Version.h"

namespace CepGen
{
  const char* version()
  {
    char* str = new char[9];
    sprintf( str, "%02u.%02u.%02u", ( cepgen_version >> 16 ) & 0xff,
                                    ( cepgen_version >>  8 ) & 0xff,
                                      cepgen_version         & 0xff );
    return str;
  }
}
