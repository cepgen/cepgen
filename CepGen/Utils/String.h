#ifndef CepGen_Utils_String_h
#define CepGen_Utils_String_h

#include <string>
#include <vector>
#include <numeric>

namespace cepgen {
  namespace utils {
    /// Format a string using a printf style format descriptor.
    std::string format(const std::string fmt, ...);
    /// Human-readable boolean printout
    std::string yesno(bool test);
    /// Boldify a string for TTY-type output streams
    /// \tparam T type of variable to be boldified
    template <typename T>
    std::string boldify(T str);
    /// TTY-type enumeration of colours
    enum class Colour {
      reset = -1,
      black = 30,
      red = 31,
      green = 32,
      yellow = 33,
      blue = 34,
      magenta = 35,
      cyan = 36,
      white = 37
    };
    enum class Modifier { reset = -1, bold = 1, dimmed = 2, italic = 3, underline = 4, blink = 5, reverse = 7 };
    /// Colourise a string for TTY-type output streams
    std::string colourise(const std::string& str, const Colour& col, const Modifier& mod = Modifier::reset);
    /// Replace all occurrences of a text by another
    size_t replace_all(std::string& str, const std::string& from, const std::string& to);
    /// Split a string according to a separation character
    std::vector<std::string> split(const std::string&, char);
    /// Merge a collection of strings in a single string
    std::string merge(const std::vector<std::string>&, const std::string&);
    /// Capitalise a string
    std::string toupper(const std::string&);
    /// Lowercase version of a string
    std::string tolower(const std::string&);
    /// Add a trailing "s" when needed
    inline const char* s(size_t num) { return (num > 1) ? "s" : ""; }
    /// Add a trailing "s" when needed
    inline std::string s(const std::string& word, size_t num, bool show_number = true) {
      return show_number ? format("%i %s%s", num, word.c_str(), (num > 1) ? "s" : "")
                         : format("%s%s", word.c_str(), (num > 1) ? "s" : "");
    }
    /// Helper to print a vector
    template <class T>
    std::string repr(const std::vector<T>& vec, const std::string& sep = ",") {
      return std::accumulate(
          std::next(vec.begin()), vec.end(), std::to_string(*vec.begin()), [&sep](std::string str, T xv) {
            return std::move(str) + sep + std::to_string(xv);
          });
    }
    /// Trim leading spaces
    void ltrim(std::string&);
    /// Trim trailing spaces
    void rtrim(std::string&);
    /// Trim leading and trailing spaces
    inline void trim(std::string& s) {
      ltrim(s);
      rtrim(s);
    }
    /// Get an environment variable
    std::string environ(const std::string&, const std::string& def = "");
  }  // namespace utils
}  // namespace cepgen

/// Provide a random number generated along a uniform distribution between 0 and 1
#define drand() (double)rand() / RAND_MAX

#endif
