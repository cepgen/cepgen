#ifndef CepGen_Utils_Filesystem_h
#define CepGen_Utils_Filesystem_h
#if __cplusplus >= 201703L
#include <filesystem>
namespace fs = std::filesystem;
#elif __cplusplus >= 201103L
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#else
#error "*** no support for filesystem! ***"
#endif

#include <string>

namespace cepgen {
  namespace utils {
    bool fileExists(const std::string& path);
  }  //namespace utils
}  // namespace cepgen

#endif
