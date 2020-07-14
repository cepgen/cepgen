#ifndef CepGen_Version_h
#define CepGen_Version_h

#include <string>

namespace cepgen
{
  struct version
  {
    /// CepGen version
    static const std::string tag;
    /// CepGen detailed version
    static const std::string extended;
  };
}

#endif
