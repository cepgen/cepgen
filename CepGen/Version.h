#ifndef CepGen_Version_h
#define CepGen_Version_h

#include <string>

namespace cepgen
{
  /// CepGen version
  /// \note Format: 0xMMmmff, with
  ///  - MM = major version
  ///  - mm = minor version
  ///  - ff = feature(s) release
  const unsigned int cepgen_version = 0x000900;
  /// Human-readable version number
  const std::string version();
}

#endif
