/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
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
  namespace utils {
    /// Convert a wide characters to a standard characters string
    std::string tostring(const std::wstring& str);
    /// Convert a wide characters to a standard characters string
    std::wstring towstring(const std::string& str);
    /// Format a string using a printf style format descriptor.
    template <typename... Args>
    std::string format(const std::string& fmt, Args... args) {
      // first check how much space is required for output buffer
      size_t size = snprintf(nullptr, 0, fmt.data(), args...) + 1;  // extra space for last '\0'
      if (size <= 0)
        return fmt;
      std::vector<char> buffer(size);
      snprintf(buffer.data(), size, fmt.data(), args...);
      return std::string(buffer.data(), buffer.data() + size - 1);  // strip last '\0'
    }
    template <typename... Args>
    std::string format(const std::wstring& fmt, Args... args) {
      return format(tostring(fmt), args...);
    }
    /// Demangle a type id if possible
    std::string demangle(const char*);
    /// Return the formatted date/time now
    std::string timeAs(const std::string& fmt);
    /// Human-readable boolean printout
    std::string yesno(bool test);
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
    std::string colourise(const std::string& str, const Colour& col, const Modifier& mod = Modifier::none);
    /// Replace all unsafe characters to build a computer-readable (and filename-safe) string
    std::string sanitise(const std::string&);
    /// Transform all emoji-like special characters into their LaTeX representation
    std::string parseSpecialChars(const std::string&);
    /// Replace all occurrences of a text by another
    size_t replace_all(std::string& str, const std::string& from, const std::string& to);
    /// Replace all occurrences of a text by another
    std::string replace_all(const std::string& str, const std::string& from, const std::string& to);
    /// Replace all occurrences of multiple texts by others
    std::string replace_all(const std::string& str, const std::vector<std::pair<std::string, std::string> >& keys);
    /// Split a string according to a separation character
    std::vector<std::string> split(const std::string&, char, bool trim = false);
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
    /// Check if a string is also an integer
    bool isInt(const std::string&);
    /// Check if a string is also a floating point number
    bool isFloat(const std::string&);
    /// Transform any type into a string
    template <typename T>
    std::string to_string(const T&);
    /// Check if a collection contains an item
    template <typename T>
    bool contains(const std::vector<T>& coll, const T& item) {
      return std::find(coll.begin(), coll.end(), item) != coll.end();
    }
    /// Check if a collection contains an item
    template <typename T>
    bool contains(const std::set<T>& coll, const T& item) {
      return std::find(coll.begin(), coll.end(), item) != coll.end();
    }
    template <typename K, typename T>
    bool contains(const std::unordered_map<K, T>& coll, const T& item) {
      return std::find_if(coll.begin(), coll.end(), [&item](const auto& kv) { return kv.second == item; }) !=
             coll.end();
    }
    /// Remove duplicates and sort a collection
    template <typename T>
    void normalise(std::vector<T>& coll);
    /// Check if all elements of a collection are uniform
    template <typename T>
    inline bool uniform(const std::vector<T>& coll) {
      return coll.size() > 1 ? coll == std::vector<T>(coll.size(), coll.at(0)) : true;
    }
    /// Capitalise a string
    std::string toupper(const std::string&);
    /// Lowercase version of a string
    std::string tolower(const std::string&);
    /// Get a (list of) substring(s) between two characters chains
    /// \param[in] beg Start delimiter of the substring(s)
    /// \param[in] end End delimiter of the substring(s)
    std::vector<std::string> between(const std::string& str, const std::string& beg, const std::string& end);
    /// Add a trailing "s" when needed
    inline const char* s(size_t num) { return (num > 1) ? "s" : ""; }
    /// Add a trailing "s" when needed
    inline std::string s(const std::string& word, float num, bool show_number = true) {
      return show_number ? format("%g %s%s", num, word.c_str(), (num > 1.) ? "s" : "")
                         : format("%s%s", word.c_str(), (num > 1.) ? "s" : "");
    }
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
    template <class T>
    inline std::string repr(const std::vector<T>& vec, const std::string& sep = ",") {
      return repr<T>(
          vec, [](const T& xv) { return std::to_string(xv); }, sep);
    }
    template <>
    inline std::string repr(const std::vector<std::string>& vec, const std::string& sep) {
      return repr<std::string>(
          vec, [](const std::string& xv) { return xv; }, sep);
    }
    std::string randomString(size_t size);
    /// Trim leading spaces
    std::string ltrim(const std::string& str);
    /// Trim trailing spaces
    std::string rtrim(const std::string& str);
    /// Trim leading and trailing spaces
    inline std::string trim(const std::string& str) { return ltrim(rtrim(str)); }
    /// Strip all special characters from string
    std::string strip(const std::string&);
    /// Check if a string starts with a given token
    bool startsWith(const std::string&, const std::string&);
    /// Describe an error code
    std::string describeError(int errnum);
  }  // namespace utils
}  // namespace cepgen

#endif
