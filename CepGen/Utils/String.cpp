#include <unistd.h>

#include <cmath>
#include <codecvt>
#include <iterator>
#include <locale>
#include <sstream>
#include <unordered_set>
#include <vector>

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  namespace utils {
    std::string yesno(bool test) { return test ? colourise("yes", Colour::green) : colourise("no", Colour::red); }

    /// String implementation of the boldification procedure
    template <>
    std::string boldify(std::string str) {
      return colourise(str, Colour::reset, Modifier::bold);
    }

    /// C-style character array implementation of the boldification procedure
    template <>
    std::string boldify(const char* str) {
      return boldify(std::string(str));
    }

    /// Unsigned long integer implementation of the boldification procedure
    template <>
    std::string boldify(unsigned long ui) {
      return boldify(std::to_string(ui));
    }

    std::string colourise(const std::string& str, const Colour& col, const Modifier& mod) {
      if (!isatty(fileno(stdout)))
        return str;
      if (mod == Modifier::reset)
        return format("\033[%dm%s\033[0m", (int)col, str.c_str());
      if (col == Colour::reset)
        return format("\033[%dm%s\033[0m", (int)mod, str.c_str());
      return format("\033[%d;%dm%s\033[0m", (int)col, (int)mod, str.c_str());
    }

    size_t replace_all(std::string& str, const std::string& from, const std::string& to) {
      size_t count = 0, pos = 0;
      while ((pos = str.find(from, pos)) != std::string::npos) {
        str.replace(pos, from.length(), to);
        pos += to.length();
        ++count;
      }
      return count;
    }

    std::vector<std::string> split(const std::string& str, char delim) {
      std::vector<std::string> out;
      std::string token;
      std::istringstream iss(str);
      while (std::getline(iss, token, delim))
        out.emplace_back(token);
      return out;
    }

    std::string merge(const std::vector<std::string>& vec, const std::string& delim) {
      if (vec.empty())
        return std::string();
      if (vec.size() == 1)
        return vec.at(0);
      std::ostringstream oss;
      std::copy(vec.begin(), std::prev(vec.end()), std::ostream_iterator<std::string>(oss, delim.c_str()));
      return oss.str() + *vec.rbegin();  // treat last one separately to drop the last delimiter
    }

    std::string toupper(const std::string& str) {
      std::string out;
      out.resize(str.size());
      std::transform(str.begin(), str.end(), out.begin(), ::toupper);
      return out;
    }

    std::string tolower(const std::string& str) {
      std::string out;
      out.resize(str.size());
      std::transform(str.begin(), str.end(), out.begin(), ::tolower);
      return out;
    }

    template <typename T>
    void normalise(std::vector<T>& coll) {
      std::unordered_set<T> set;
      for (const auto& it : coll)
        set.insert(it);
      coll.assign(set.begin(), set.end());
      std::sort(coll.begin(), coll.end());
    }
    template void normalise(std::vector<std::string>&);

    void ltrim(std::string& s) {
      s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) { return !std::isspace(ch); }));
    }

    void rtrim(std::string& s) {
      s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) { return !std::isspace(ch); }).base(), s.end());
    }

    std::string strip(const std::string& str) {
      auto out = str;
      out.resize(std::remove_if(out.begin(), out.end(), [](char x) { return !std::isalnum(x) && !std::isspace(x); }) -
                 out.begin());
      return out;
    }

    std::string tostring(const std::wstring& str) {
      typedef std::codecvt_utf8_utf16<wchar_t> convert_type;
      std::wstring_convert<convert_type, wchar_t> converter;
      return converter.to_bytes(str);
    }
    std::wstring towstring(const std::string& str) {
      typedef std::codecvt_utf8_utf16<wchar_t> convert_type;
      std::wstring_convert<convert_type, wchar_t> converter;
      return converter.from_bytes(str);
    }

    namespace env {
      std::string get(const std::string& var, const std::string& def) {
        const auto out = std::getenv(var.c_str());
        if (!out)
          return def;
        return std::string(out);
      }

      void set(const std::string& var, const std::string& value) { setenv(var.c_str(), value.c_str(), 1); }

#ifdef _WIN32
      static constexpr const char* PATH_DELIM = ";";
#else
      static constexpr const char* PATH_DELIM = ":";
#endif

      void append(const std::string& var, const std::string& value) {
        auto env = split(get(var), PATH_DELIM[0]);
        env.emplace_back(value);
        normalise(env);
        setenv(var.c_str(), merge(env, PATH_DELIM).c_str(), 1);
      }

      void unset(const std::string& var) { unsetenv(var.c_str()); }
    }  // namespace env
  }    // namespace utils
}  // namespace cepgen
