#include "CepGen/Version.h"
#include "CepGen/Core/utils.h"

namespace cepgen
{
  const std::string version()
  {
    return Form( "%02u.%02u.%02u", ( cepgen_version >> 16 ) & 0xff,
                                   ( cepgen_version >>  8 ) & 0xff,
                                     cepgen_version         & 0xff );
  }
}
