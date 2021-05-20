#include "CepGen/Utils/Filesystem.h"
#include <fstream>

namespace cepgen {
  namespace utils {
    bool fileExists(const std::string& path) { return std::ifstream(path).good(); }
  }  // namespace utils
}  // namespace cepgen
