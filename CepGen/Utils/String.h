/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2024  Laurent Forthomme
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CepGen_Utils_String_h
#define CepGen_Utils_String_h

#include <algorithm>
#include <functional>
#include <numeric>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

namespace cepgen {
  class Limits;
  class ParametersList;
  namespace utils {
    /// Transform any type into a string
    template <typename T>
    inline std::string toString(const T& obj) {
      return std::to_string(obj);
    }
    /// Specialisation of string conversion of limits
    template <>
    std::string toString(const Limits&);
    /// Specialisation of string conversion of parameters list
    template <>
    std::string toString(const ParametersList&);
    /// Trivial specialisation of string conversion of string
    template <>
    inline std::string toString(const std::string& obj) {
      return obj;
    }
    template <>
    std::string toString(const std::wstring&);       ///< Convert a wide characters to a standard characters string
    std::wstring toWstring(const std::string& str);  ///< Convert a standard characters to a wide characters string
    std::string toCamelCase(const std::string&, bool lower = true);  ///< Convert any case into a camelCase string
    bool isInt(const std::string&);                                  ///< Check if a string is also an integer
    bool isFloat(const std::string&);  ///< Check if a string is also a floating point number

    /// Format a string using a printf style format descriptor.
    template <typename... Args>
    inline std::string format(const std::string& fmt, Args... args) {
      // first check how much space is required for output buffer
      size_t size = snprintf(nullptr, 0, fmt.data(), args...) + 1;  // extra space for last '\0'
      if (size <= 0)
        return fmt;
      std::vector<char> buffer(size);
      snprintf(buffer.data(), size, fmt.data(), args...);
      return std::string(buffer.data(), buffer.data() + size - 1);  // strip last '\0'
    }
    /// Format a wide string using a printf style format descriptor
    template <typename... Args>
    inline std::string format(const std::wstring& fmt, Args... args) {
      return format(toString(fmt), args...);
    }
    std::string demangle(const char*);           ///< Demangle a type id if possible
    std::string timeAs(const std::string& fmt);  ///< Return the formatted date/time now
    std::string yesno(bool test);                ///< Human-readable boolean printout
    /// Boldify a string for TTY-type output streams
    /// \tparam T type of variable to be boldified
    template <typename T>
    std::string boldify(T str);
    /// TTY-type enumeration of colours
    enum class Colour {
      none = -1,
      reset = 0,
      black = 30,
      red = 31,
      green = 32,
      yellow = 33,
      blue = 34,
      magenta = 35,
      cyan = 36,
      white = 37
    };
    std::ostream& operator<<(std::ostream&, const Colour&);
    enum struct Modifier : int16_t {
      none = -1,
      reset = 0,
      bold = 1,
      dimmed = 1 << 1,
      italic = 1 << 2,
      underline = 1 << 3,
      blink = 1 << 4,
      reverse = 1 << 6
    };
    std::ostream& operator<<(std::ostream&, const Modifier&);
    Modifier operator|(const Modifier&, const Modifier&);
    /// Colourise a string for TTY-type output streams
    std::string colourise(const std::string&, const Colour&, const Modifier& = Modifier::none);
    /// Replace all unsafe characters to build a computer-readable (and filename-safe) string
    std::string sanitise(const std::string&);
    /// Transform all emoji-like special characters into their LaTeX representation
    std::string parseSpecialChars(const std::string&);
    /// Replace all occurrences of a text by another
    size_t replaceAll(std::string& str, const std::string& from, const std::string& to);
    /// Replace all occurrences of a text by another
    std::string replaceAll(const std::string& str, const std::string& from, const std::string& to);
    /// Replace all occurrences of multiple texts by others
    std::string replaceAll(const std::string& str, const std::vector<std::pair<std::string, std::string> >& keys);
    /// Split a string according to a separation character
    std::vector<std::string> split(const std::string&, char, bool trim = false);
    /// Merge a a printable type in a single string
    template <typename T>
    inline std::string merge(const T& val, const std::string&) {
      return toString(val);
    }
    /// Merge a collection of a printable type in a single string
    template <typename T>
    std::string merge(const std::vector<T>&, const std::string&);
    /// Merge a collection of collections of a printable type in a single string
    template <typename T>
    std::string merge(const std::vector<std::vector<T> >&, const std::string&);
    /// Merge a collection of a printable type in a single string
    template <typename T, size_t N>
    inline std::string merge(const std::array<T, N>& arr, const std::string& delim) {
      return merge(std::vector<T>(arr.begin(), arr.end()), delim);
    }
    /// Trivial dimension-1 "merger" for string input
    inline std::string merge(const std::string& val, const std::string&) { return val; }
    /// Trivial dimension-1 "merger" for parameters list input
    std::string merge(const ParametersList&, const std::string&);
    /// Trivial dimension-1 "merger" for limits input
    std::string merge(const Limits&, const std::string&);

    std::string toUpper(const std::string&);  ///< Capitalise a string
    std::string toLower(const std::string&);  ///< Lowercase version of a string

    /// Get a (list of) substring(s) between two characters chains
    /// \param[in] beg Start delimiter of the substring(s)
    /// \param[in] end End delimiter of the substring(s)
    std::vector<std::string> between(const std::string& str, const std::string& beg, const std::string& end);
    std::string s(const std::string&, float, bool = true);  ///< Add a trailing "s" when needed

    /// Helper to print a vector
    template <class T>
    inline std::string repr(const std::vector<T>& vec,
                            const std::function<std::string(const T&)>& printer,
                            const std::string& sep = ",") {
      if (vec.empty())
        return "{}";
      return std::accumulate(
          std::next(vec.begin()), vec.end(), printer(*vec.begin()), [&printer, &sep](std::string str, const T& xv) {
            return std::move(str) + sep + printer(xv);
          });
    }
    /// Helper to print a vector
    template <class T>
    inline std::string repr(const std::vector<T>& vec, const std::string& sep = ",") {
      return repr<T>(
          vec, [](const T& xv) { return toString(xv); }, sep);
    }
    std::string randomString(size_t size);      ///< Generate a random string of a given size
    std::string ltrim(const std::string& str);  ///< Trim leading spaces
    std::string rtrim(const std::string& str);  ///< Trim trailing spaces
    inline std::string trim(const std::string& str) { return ltrim(rtrim(str)); }  ///< Trim leading and trailing spaces
    std::string strip(const std::string&);                    ///< Strip all special characters from string
    bool startsWith(const std::string&, const std::string&);  ///< Check if a string starts with a given token
    bool endsWith(const std::string&, const std::string&);    ///< Check if a string ends with a given token
    std::string describeError(int errnum);                    ///< Describe an error code

  }  // namespace utils
}  // namespace cepgen

#endif
