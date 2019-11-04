#include "CepGen/Version.h"
#include "CepGen/Utils/String.h"

namespace cepgen
{
  const std::string version()
  {
    return Form( "%u.%u.%u", ( cepgen_version >> 16 ) & 0xff,
                             ( cepgen_version >>  8 ) & 0xff,
                               cepgen_version         & 0xff );
  }
}
