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
    /// Check if the file exists
    inline bool fileExists(const std::string& path) { return fs::exists(path); }
    /// Small utility to retrieve the extension of a filename
    inline std::string fileExtension(const std::string& file) { return fs::path(file).extension(); }
  }  //namespace utils
}  // namespace cepgen

#endif
