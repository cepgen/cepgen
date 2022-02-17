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

#include <unistd.h>

#include <cmath>
#include <codecvt>
#include <locale>
#include <unordered_set>
#include <vector>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  namespace utils {
    std::string yesno(bool test) { return test ? colourise("true", Colour::green) : colourise("false", Colour::red); }

    /// String implementation of the boldification procedure
    template <>
    std::string boldify(std::string str) {
      return colourise(str, Colour::none, Modifier::bold);
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

    Modifier operator|(const Modifier& lhs, const Modifier& rhs) {
      std::bitset<7> mod1((int)lhs), mod2((int)rhs);
      return (Modifier)(mod1 | mod2).to_ulong();
    }

    std::string colourise(const std::string& str, const Colour& col, const Modifier& mod) {
      if (!isatty(fileno(stdout)))
        return str;
      std::string out;
      auto get_mod_str = [](const Colour& col, const Modifier& mod) -> std::string {
        std::string mod_str("\033[");
        if (col != Colour::none)
          mod_str += std::to_string((int)col);
        if (mod > Modifier::reset)
          for (size_t i = 0; i < 7; ++i)
            if (((uint16_t)mod >> i) & 0x1)
              mod_str += ";" + std::to_string(i + 1);
        return mod_str + "m";
      };
      out = get_mod_str(col, mod);
      out += str;
      out += get_mod_str(Colour::reset, Modifier::reset);
      return out;
    }

    std::string parseSpecialChars(const std::string& str) {
      return replace_all(
          str, {{"Α", "\\Alpha"},      {"Β", "\\Beta"},      {"Χ", "\\Chi"},     {"Δ", "\\Delta"},   {"Ε", "\\Epsilon"},
                {"Φ", "\\Phi"},        {"Γ", "\\Gamma"},     {"Η", "\\Eta"},     {"Ι", "\\Iota"},    {"Κ", "\\Kappa"},
                {"Λ", "\\Lambda"},     {"Μ", "\\Mu"},        {"Ν", "\\Nu"},      {"Ο", "\\Omicron"}, {"Π", "\\Pi"},
                {"Θ", "\\Theta"},      {"Ρ", "\\Rho"},       {"Σ", "\\Sigma"},   {"Τ", "\\Tau"},     {"Υ", "\\Upsilon"},
                {"Ω", "\\Omega"},      {"Ξ", "\\Xi"},        {"Ψ", "\\Psi"},     {"Ζ", "\\Zeta"},    {"α", "\\alpha"},
                {"β", "\\beta"},       {"χ", "\\Chi"},       {"δ", "\\delta"},   {"ε", "\\epsilon"}, {"ɸ", "\\phi"},
                {"γ", "\\gamma"},      {"η", "\\eta"},       {"ι", "\\iota"},    {"κ", "\\kappa"},   {"λ", "\\lambda"},
                {"μ", "\\mu"},         {"ν", "\\nu"},        {"ο", "\\omicron"}, {"π", "\\pi"},      {"θ", "\\theta"},
                {"ρ", "\\rho"},        {"σ", "\\sigma"},     {"τ", "\\tau"},     {"υ", "\\upsilon"}, {"ω", "\\omega"},
                {"ξ", "\\xi"},         {"ψ", "\\psi"},       {"ζ", "\\zeta"},    {"⁺", "^{+}"},      {"¯", "^{-}"},
                {"→", "\\rightarrow"}, {"←", "\\leftarrow"}, {"↝ ", "\\leadsto"}});
    }

    std::string sanitise(const std::string& str) {
      return replace_all(str, {{")", ""}, {"(", "_"}, {"{", "_"}, {",", "_"}, {":", "_"}});
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

    std::string replace_all(const std::string& str, const std::string& from, const std::string& to) {
      auto out{str};
      if (replace_all(out, from, to) == 0)
        CG_DEBUG("replace_all") << "No occurences of {" << from << "} found in input string.";
      return out;
    }

    std::string replace_all(const std::string& str, const std::vector<std::pair<std::string, std::string> >& keys) {
      auto out{str};
      for (const auto& key : keys)
        replace_all(out, key.first, key.second);
      CG_DEBUG_LOOP("replace_all").log([&keys, &out](auto& log) {
        log << "Values to be replaced: ";
        for (const auto& key : keys)
          log << "\n\t{\"" << key.first << "\" -> \"" << key.second << "\"}";
        log << "\n-> output: \"" << out << "\".";
      });
      return out;
    }

    std::vector<std::string> split(const std::string& str, char delim) {
      std::vector<std::string> out;
      if (str.empty())
        return out;
      std::string token;
      std::istringstream iss(str);
      while (std::getline(iss, token, delim))
        out.emplace_back(token);
      return out;
    }

    template <typename T>
    std::string merge(const std::vector<T>& vec, const std::string& delim) {
      if (vec.empty())
        return std::string();
      std::ostringstream oss;
      std::for_each(vec.begin(), vec.end(), [&oss, &delim, sep = std::string()](const auto& val) mutable {
        oss << sep << val;
        sep = delim;
      });
      return oss.str();
    }

    template std::string merge<std::string>(const std::vector<std::string>&, const std::string&);
    template std::string merge<int>(const std::vector<int>&, const std::string&);
    template std::string merge<double>(const std::vector<double>&, const std::string&);
    template std::string merge<ParametersList>(const std::vector<ParametersList>&, const std::string&);

    bool isNumber(const std::string& str) {
      return !str.empty() &&
             std::find_if(str.begin(), str.end(), [](unsigned char c) { return !std::isdigit(c); }) == str.end();
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

    std::string ltrim(const std::string& str) {
      auto out{str};
      out.erase(out.begin(), std::find_if(out.begin(), out.end(), [](unsigned char ch) { return !std::isspace(ch); }));
      return out;
    }

    std::string rtrim(const std::string& str) {
      auto out{str};
      out.erase(std::find_if(out.rbegin(), out.rend(), [](unsigned char ch) { return !std::isspace(ch); }).base(),
                out.end());
      return out;
    }

    std::string strip(const std::string& str) {
      auto out{str};
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
      static constexpr char PATH_DELIM = ';';
#else
      static constexpr char PATH_DELIM = ':';
#endif

      void append(const std::string& var, const std::string& value) {
        auto env = split(get(var), PATH_DELIM);
        env.emplace_back(value);
        normalise(env);
        setenv(var.c_str(), merge(env, std::string(1, PATH_DELIM)).c_str(), 1);
      }

      void unset(const std::string& var) { unsetenv(var.c_str()); }
    }  // namespace env

    std::string describeError(int errnum) {
      std::array<char, 50> buff;
      return std::to_string(errnum) + " (" + std::string(strerror_r(errnum, buff.data(), buff.size())) + ")";
    }

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define DEF_COLOUR(col) \
  case Colour::col:     \
    return os << colourise(TOSTRING(col), Colour::col);
#define DEF_MODIFIER(mod) \
  case Modifier::mod:     \
    return os << colourise(TOSTRING(mod), Colour::none, Modifier::mod);

    std::ostream& operator<<(std::ostream& os, const Colour& col) {
      switch (col) {
        DEF_COLOUR(reset)
        DEF_COLOUR(black)
        DEF_COLOUR(red)
        DEF_COLOUR(green)
        DEF_COLOUR(yellow)
        DEF_COLOUR(blue)
        DEF_COLOUR(magenta)
        DEF_COLOUR(cyan)
        DEF_COLOUR(white)
        default:
          DEF_COLOUR(none)
      }
    }

    std::ostream& operator<<(std::ostream& os, const Modifier& mod) {
      switch (mod) {
        DEF_MODIFIER(reset)
        DEF_MODIFIER(bold)
        DEF_MODIFIER(dimmed)
        DEF_MODIFIER(italic)
        DEF_MODIFIER(underline)
        DEF_MODIFIER(blink)
        DEF_MODIFIER(reverse)
        default:
          DEF_MODIFIER(none)
      }
    }

#undef DEF_COLOUR
#undef DEF_MODIFIER
  }  // namespace utils
}  // namespace cepgen
