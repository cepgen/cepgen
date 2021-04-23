#ifndef CepGen_Version_h
#define CepGen_Version_h

#include <string>

namespace cepgen {
  /// Collection of CepGen version information handlers
  struct version {
    /// CepGen version
    static const std::string tag;
    /// CepGen detailed version
    static const std::string extended;
    /// CepGen banner
    static const std::string banner;
  };
}  // namespace cepgen

#endif
