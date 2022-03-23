/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2018-2021  Laurent Forthomme
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

#include <climits>
#include <iomanip>
#include <regex>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/String.h"

#define IMPL_TYPE_GET(type, coll, name)                                                                  \
  template <>                                                                                            \
  type ParametersList::get<type>(const std::string& key, const type& def) const {                        \
    if (coll.count(key) > 0)                                                                             \
      return coll.at(key);                                                                               \
    CG_DEBUG("ParametersList") << "Failed to retrieve " << name << " parameter with key=" << key << ". " \
                               << "Default value: " << def << ".";                                       \
    return def;                                                                                          \
  }

#define IMPL_TYPE_SET(type, coll, name)                                                                             \
  template <>                                                                                                       \
  bool ParametersList::has<type>(const std::string& key) const {                                                    \
    return coll.count(key) != 0;                                                                                    \
  }                                                                                                                 \
  template <>                                                                                                       \
  ParametersList& ParametersList::set<type>(const std::string& key, const type& value) {                            \
    coll[key] = value;                                                                                              \
    return *this;                                                                                                   \
  }                                                                                                                 \
  template <>                                                                                                       \
  type& ParametersList::operator[]<type>(const std::string& key) {                                                  \
    return coll[key];                                                                                               \
  }                                                                                                                 \
  template <>                                                                                                       \
  std::vector<std::string> ParametersList::keysOf<type>() const {                                                   \
    std::vector<std::string> out;                                                                                   \
    std::transform(coll.begin(), coll.end(), std::back_inserter(out), [](const auto& pair) { return pair.first; }); \
    return out;                                                                                                     \
  }

#define IMPL_TYPE_ALL(type, coll, name) \
  IMPL_TYPE_GET(type, coll, #name)      \
  IMPL_TYPE_SET(type, coll, #name)

namespace cepgen {
  const std::string ParametersList::MODULE_NAME = "mod_name";
  const std::regex kFloatRegex("[+-]?([0-9]+)[.EeDd][+-]?([0-9]*)?|[.][0-9]+", std::regex_constants::extended);

  ParametersList::ParametersList(const ParametersList& oth)
      : param_values_(oth.param_values_),
        bool_values_(oth.bool_values_),
        int_values_(oth.int_values_),
        ulong_values_(oth.ulong_values_),
        dbl_values_(oth.dbl_values_),
        str_values_(oth.str_values_),
        lim_values_(oth.lim_values_),
        vec_int_values_(oth.vec_int_values_),
        vec_dbl_values_(oth.vec_dbl_values_),
        vec_str_values_(oth.vec_str_values_),
        vec_param_values_(oth.vec_param_values_) {}

  bool ParametersList::operator==(const ParametersList& oth) const {
    // only ensure the keys are identical
    /// \note this might be a bit too loose for more advanced usages
    if (keys() != oth.keys())
      return false;
    return true;
  }

  ParametersList& ParametersList::operator+=(const ParametersList& oth) {
    // ensure the two collections are not identical
    if (*this == oth)
      return *this;
    // then check if any key of the other collection is lready present in the list
    std::vector<std::string> keys_erased;
    for (const auto& key : oth.keys()) {
      if (has<ParametersList>(key)) {
        // do not remove a duplicate parameters collection if they are not strictly identical ;
        // will concatenate its values with the other object's
        if (get<ParametersList>(key) == oth.get<ParametersList>(key) && erase(key) > 0)
          keys_erased.emplace_back(key);
      } else if (erase(key) > 0)
        // any other duplicate key is just replaced
        keys_erased.emplace_back(key);
    }
    if (!keys_erased.empty())
      CG_DEBUG("ParametersList") << utils::s("key", keys_erased.size(), true) << " erased: " << keys_erased << ".";
    //--- concatenate all typed lists
    bool_values_.insert(oth.bool_values_.begin(), oth.bool_values_.end());
    int_values_.insert(oth.int_values_.begin(), oth.int_values_.end());
    vec_int_values_.insert(oth.vec_int_values_.begin(), oth.vec_int_values_.end());
    dbl_values_.insert(oth.dbl_values_.begin(), oth.dbl_values_.end());
    vec_dbl_values_.insert(oth.vec_dbl_values_.begin(), oth.vec_dbl_values_.end());
    str_values_.insert(oth.str_values_.begin(), oth.str_values_.end());
    vec_str_values_.insert(oth.vec_str_values_.begin(), oth.vec_str_values_.end());
    vec_param_values_.insert(oth.vec_param_values_.begin(), oth.vec_param_values_.end());
    lim_values_.insert(oth.lim_values_.begin(), oth.lim_values_.end());
    // special case for parameters collection: concatenate values instead of full containers
    for (const auto& par : oth.param_values_)
      // if the two parameters list are modules, and do not have the same name,
      // simply replace the old one with the new parameters list
      if (param_values_[par.first].getString(ParametersList::MODULE_NAME) ==
          par.second.getString(ParametersList::MODULE_NAME))
        param_values_[par.first] += par.second;
      else
        param_values_[par.first] = par.second;
    return *this;
  }

  ParametersList ParametersList::operator+(const ParametersList& oth) const {
    ParametersList out = *this;
    out += oth;
    return out;
  }

  ParametersList& ParametersList::feed(const std::string& raw_args) {
    auto raw_list = raw_args;
    const auto raw_list_stripped = utils::between(raw_list, "{", "}");
    if (raw_list_stripped.size() == 1 && raw_list == "{" + raw_list_stripped[0] + "}")
      raw_list = raw_list_stripped[0];
    // first preprocess the arguments list to isolate all comma-separated arguments
    std::vector<std::string> list;
    std::vector<std::string> buf;
    short num_open_braces = 0;
    for (const auto& item : utils::split(raw_list, ',')) {
      buf.emplace_back(item);
      num_open_braces += std::count(item.begin(), item.end(), '{') - std::count(item.begin(), item.end(), '}');
      if (num_open_braces <= 0) {
        list.emplace_back(utils::merge(buf, ","));
        buf.clear();
      }
    }
    CG_DEBUG("ParametersList:feed") << "Parsed arguments: " << list << ", raw list: " << raw_list
                                    << " (splitted: " << utils::split(raw_list, ',') << "), "
                                    << "{-} imbalance: " << num_open_braces << ".";
    if (num_open_braces != 0)
      throw CG_ERROR("ParametersList:feed") << "Invalid string to be parsed as a parameters list!\n\t"
                                            << "Open-closed braces imbalance: " << num_open_braces << "\n\t"
                                            << "Raw list: " << raw_list << "\n\t"
                                            << "Resulting list: " << list << ", buffer: " << buf << ".";
    // now loop through all unpacked arguments
    for (const auto& arg : list) {
      // browse through the parameters hierarchy
      auto cmd = utils::split(arg, '/');
      if (cmd.size() > 1) {  // sub-parameters word found
        operator[]<ParametersList>(cmd.at(0)).feed(
            utils::merge(std::vector<std::string>(cmd.begin() + 1, cmd.end()), "/"));
        continue;
      }

      // from this moment on, a "key=value" or "key(=true)" was found
      const auto& subplist = utils::between(arg, "{", "}");
      if (!subplist.empty()) {
        for (const auto& subp : subplist)
          feed(subp);
        return *this;
      }
      const auto& word = cmd.at(0);
      auto words = utils::split(arg, '=');
      auto key = words.at(0);
      if (erase(key) > 0)
        CG_DEBUG("ParametersList:feed") << "Replacing key='" << key << "' with a new value.";
      if (key == "name")
        key = ParametersList::MODULE_NAME;
      if (words.size() == 1)  // basic key=true
        set<bool>(key, true);
      else if (words.size() == 2) {  // basic key=value
        const auto value = words.at(1);
        try {
          if (std::regex_match(value, kFloatRegex))
            set<double>(key, std::stod(value));
          else
            set<int>(key, std::stoi(value));
        } catch (const std::invalid_argument&) {
          const auto value_lc = utils::tolower(value);
          if (value_lc == "off" || value_lc == "no" || value_lc == "false")
            set<bool>(key, false);
          else if (value_lc == "on" || value_lc == "yes" || value_lc == "true")
            set<bool>(key, true);
          else
            set<std::string>(key, value);
        }
      } else
        throw CG_FATAL("ParametersList:feed") << "Invalid key=value unpacking: " << word << "!";
    }
    return *this;
  }

  size_t ParametersList::erase(const std::string& key) {
    size_t out = 0ull;
    if (bool_values_.count(key) > 0)
      out += bool_values_.erase(key);
    if (int_values_.count(key) > 0)
      out += int_values_.erase(key);
    if (ulong_values_.count(key) > 0)
      out += ulong_values_.erase(key);
    if (dbl_values_.count(key) > 0)
      out += dbl_values_.erase(key);
    if (str_values_.count(key) > 0)
      out += str_values_.erase(key);
    if (lim_values_.count(key) > 0)
      out += lim_values_.erase(key);
    if (param_values_.count(key) > 0)
      out += param_values_.erase(key);
    if (vec_int_values_.count(key) > 0)
      out += vec_int_values_.erase(key);
    if (vec_dbl_values_.count(key) > 0)
      out += vec_dbl_values_.erase(key);
    if (vec_str_values_.count(key) > 0)
      out += vec_str_values_.erase(key);
    if (vec_param_values_.count(key) > 0)
      out += vec_param_values_.erase(key);
    return out;
  }

  bool ParametersList::empty() const { return keys(false).empty(); }

  std::ostream& operator<<(std::ostream& os, const ParametersList& params) {
    params.print(os);
    return os;
  }

  const ParametersList& ParametersList::print(std::ostream& os) const {
    if (empty()) {
      os << "{}";
      return *this;
    }
    std::string sep;
    const auto& plist_name = getString(ParametersList::MODULE_NAME);
    if (!plist_name.empty()) {
      auto mod_name = has<std::string>(MODULE_NAME) ? "\"" + plist_name + "\"" : plist_name;
      os << "Module(" << mod_name, sep = ", ";
    } else
      os << "Parameters(";
    for (const auto& key : keys(false))
      os << sep << key << "=" << getString(key, true), sep = ", ";
    os << ")";
    return *this;
  }

  std::vector<std::string> ParametersList::keys(bool name_key) const {
    std::vector<std::string> out{};
    auto key = [](const auto& p) { return p.first; };
    std::transform(bool_values_.begin(), bool_values_.end(), std::back_inserter(out), key);
    std::transform(int_values_.begin(), int_values_.end(), std::back_inserter(out), key);
    std::transform(ulong_values_.begin(), ulong_values_.end(), std::back_inserter(out), key);
    std::transform(vec_int_values_.begin(), vec_int_values_.end(), std::back_inserter(out), key);
    std::transform(dbl_values_.begin(), dbl_values_.end(), std::back_inserter(out), key);
    std::transform(param_values_.begin(), param_values_.end(), std::back_inserter(out), key);
    std::transform(vec_dbl_values_.begin(), vec_dbl_values_.end(), std::back_inserter(out), key);
    std::transform(str_values_.begin(), str_values_.end(), std::back_inserter(out), key);
    std::transform(lim_values_.begin(), lim_values_.end(), std::back_inserter(out), key);
    std::transform(vec_str_values_.begin(), vec_str_values_.end(), std::back_inserter(out), key);
    std::transform(vec_param_values_.begin(), vec_param_values_.end(), std::back_inserter(out), key);
    if (!name_key) {
      const auto it_name = std::find(out.begin(), out.end(), MODULE_NAME);
      if (it_name != out.end())
        out.erase(it_name);
    }
    std::sort(out.begin(), out.end());
    return out;
  }

  std::string ParametersList::getString(const std::string& key, bool wrap) const {
    auto wrap_val = [&wrap](const auto& val, const std::string& type) -> std::string {
      std::ostringstream os;
      if (type == "float" || type == "vfloat")
        os << std::fixed;
      os << val;
      return (wrap ? type + "(" : "") + (type == "bool" ? utils::yesno(std::stoi(os.str())) : os.str()) +
             (wrap ? ")" : "");
    };
    auto wrap_coll = [&wrap, &wrap_val](const auto& coll, const std::string& type) -> std::string {
      return wrap_val(utils::merge(coll, ", "), type);
    };
    std::ostringstream os;
    if (has<ParametersList>(key))
      os << get<ParametersList>(key);
    else if (has<bool>(key))
      os << wrap_val(get<bool>(key), "bool");
    else if (has<int>(key))
      os << wrap_val(get<int>(key), "int");
    else if (has<unsigned long long>(key))
      os << wrap_val(get<unsigned long long>(key), "ulong");
    else if (has<double>(key))
      os << wrap_val(get<double>(key), "float");
    else if (has<std::string>(key))
      os << wrap_val(get<std::string>(key), "str");
    else if (has<Limits>(key))
      os << wrap_val(get<Limits>(key), "Limits");
    else if (has<std::vector<ParametersList> >(key))
      os << wrap_coll(get<std::vector<ParametersList> >(key), "VParams");
    else if (has<std::vector<int> >(key))
      os << wrap_coll(get<std::vector<int> >(key), "vint");
    else if (!has<Limits>(key) && has<std::vector<double> >(key))
      os << wrap_coll(get<std::vector<double> >(key), "vfloat");
    else if (has<std::vector<std::string> >(key))
      os << wrap_coll(get<std::vector<std::string> >(key), "vstr");
    return os.str();
  }

  ParametersList& ParametersList::rename(const std::string& old_key, const std::string& new_key) {
    if (has<bool>(old_key))
      set(new_key, get<bool>(old_key)).erase(old_key);
    if (has<int>(old_key))
      set(new_key, get<int>(old_key)).erase(old_key);
    if (has<unsigned long long>(old_key))
      set(new_key, get<unsigned long long>(old_key)).erase(old_key);
    if (has<double>(old_key))
      set(new_key, get<double>(old_key)).erase(old_key);
    if (has<std::string>(old_key))
      set(new_key, get<std::string>(old_key)).erase(old_key);
    if (has<Limits>(old_key))
      set(new_key, get<Limits>(old_key)).erase(old_key);
    if (has<ParametersList>(old_key))
      set(new_key, get<ParametersList>(old_key)).erase(old_key);
    if (has<std::vector<int> >(old_key))
      set(new_key, get<std::vector<int> >(old_key)).erase(old_key);
    if (has<std::vector<double> >(old_key))
      set(new_key, get<std::vector<double> >(old_key)).erase(old_key);
    if (has<std::vector<std::string> >(old_key))
      set(new_key, get<std::vector<std::string> >(old_key)).erase(old_key);
    if (has<std::vector<ParametersList> >(old_key))
      set(new_key, get<std::vector<ParametersList> >(old_key)).erase(old_key);
    return *this;
  }

  std::string ParametersList::serialise() const {
    std::ostringstream out;
    std::string sep;
    for (const auto& key : keys(true)) {
      out << sep << key;
      if (has<ParametersList>(key)) {
        const auto& plist = get<ParametersList>(key);
        out << "/";
        if (plist.keys().size() > 1)
          out << "{";
        out << plist.serialise();
        if (plist.keys().size() > 1)
          out << "}";
      } else
        out << "=" << getString(key, false);
      sep = ",";
    }
    return out.str();
  }

  //------------------------------------------------------------------
  // default template (placeholders)
  //------------------------------------------------------------------

  template <typename T>
  bool ParametersList::has(const std::string& key) const {
    throw CG_FATAL("ParametersList") << "Invalid type for key=" << key << "!";
  }

  template <typename T>
  T ParametersList::get(const std::string& key, const T& def) const {
    throw CG_FATAL("ParametersList") << "Invalid type retrieved for key=" << key << "!";
  }

  template <typename T>
  T& ParametersList::operator[](const std::string& key) {
    throw CG_FATAL("ParametersList") << "Invalid type retrieved for key=" << key << "!";
  }

  template <typename T>
  ParametersList& ParametersList::set(const std::string& key, const T& value) {
    throw CG_FATAL("ParametersList") << "Invalid type to be set for key=" << key << "!";
  }

  //------------------------------------------------------------------
  // sub-parameters-type attributes
  //------------------------------------------------------------------

  IMPL_TYPE_ALL(ParametersList, param_values_, "parameters")
  IMPL_TYPE_ALL(bool, bool_values_, "boolean")

  IMPL_TYPE_SET(int, int_values_, "integer")
  template <>
  int ParametersList::get<int>(const std::string& key, const int& def) const {
    if (has<int>(key))
      return int_values_.at(key);
    if (has<unsigned long long>(key)) {
      const auto ulong_val = ulong_values_.at(key);
      if (ulong_val >= INT_MAX)
        CG_WARNING("ParametersList:get")
            << "Trying to retrieve a (too) long unsigned integer with an integer getter. Please fix your code.";
      return (int)ulong_val;
    }
    return def;
  }

  IMPL_TYPE_SET(unsigned long long, ulong_values_, "unsigned long integer")
  template <>
  unsigned long long ParametersList::get<unsigned long long>(const std::string& key,
                                                             const unsigned long long& def) const {
    if (has<unsigned long long>(key))
      return ulong_values_.at(key);
    if (has<int>(key)) {
      const auto& int_val = int_values_.at(key);
      if (int_val < 0)
        CG_WARNING("ParametersList:get")
            << "Trying to retrieve a negative-value integer with an unsigned long getter. Please fix your code.";
      return int_val;
    }
    return def;
  }

  IMPL_TYPE_ALL(double, dbl_values_, "floating number")
  IMPL_TYPE_ALL(std::string, str_values_, "string")
  IMPL_TYPE_ALL(std::vector<int>, vec_int_values_, "vector of integers")
  IMPL_TYPE_ALL(std::vector<double>, vec_dbl_values_, "vector of floating numbers")
  IMPL_TYPE_ALL(std::vector<std::string>, vec_str_values_, "vector of strings")
  IMPL_TYPE_ALL(std::vector<ParametersList>, vec_param_values_, "vector of parameters")

  //------------------------------------------------------------------
  // limits-type attributes
  //------------------------------------------------------------------

  template <>
  bool ParametersList::has<Limits>(const std::string& key) const {
    if (lim_values_.count(key) != 0)
      return true;
    if (dbl_values_.count(key + "min") || dbl_values_.count(key + "max"))
      return true;
    return false;
  }

  template <>
  ParametersList& ParametersList::set<Limits>(const std::string& key, const Limits& value) {
    if (vec_dbl_values_.count(key))
      vec_dbl_values_.erase(key);
    lim_values_[key] = value;
    return *this;
  }

  template <>
  inline Limits& ParametersList::operator[]<Limits>(const std::string& key) {
    return lim_values_[key];
  }

  template <>
  Limits ParametersList::get<Limits>(const std::string& key, const Limits& def) const {
    // first try to find Limits object in collections
    auto val = std::find_if(lim_values_.begin(), lim_values_.end(), [&key](const auto& kv) { return kv.first == key; });
    if (val != lim_values_.end())
      return val->second;
    // Limits object not found ; still trying to build it from (min/max) attributes
    Limits buf;
    fill<double>(key + "min", buf.min());
    fill<double>(key + "max", buf.max());
    if (buf.valid())
      return buf.validate();
    // nothing found ; returning default
    CG_DEBUG("ParametersList") << "Failed to retrieve limits parameter with key=" << key << ". "
                               << "Default value: " << def << ".";
    return def;
  }

  template <>
  const ParametersList& ParametersList::fill<Limits>(const std::string& key, Limits& value) const {
    if (has<double>(key + "min") || has<double>(key + "max")) {
      fill<double>(key + "min", value.min());
      fill<double>(key + "max", value.max());
      return *this;
    }
    if (has<Limits>(key)) {
      const auto& lim = get<Limits>(key);
      if (lim.hasMin())
        value.min() = lim.min();
      if (lim.hasMax())
        value.max() = lim.max();
      return *this;
    }
    return *this;
  }

  template <>
  std::vector<std::string> ParametersList::keysOf<Limits>() const {
    std::vector<std::string> out;
    std::transform(
        lim_values_.begin(), lim_values_.end(), std::back_inserter(out), [](const auto& pair) { return pair.first; });
    return out;
  }

  //------------------------------------------------------------------
  // particle properties-type attributes
  //------------------------------------------------------------------

  /// Check if a particle properties object is handled
  template <>
  bool ParametersList::has<ParticleProperties>(const std::string& key) const {
    return has<ParametersList>(key) || has<int>(key);
  }

  /// Get a particle properties object
  template <>
  ParticleProperties ParametersList::get<ParticleProperties>(const std::string& key,
                                                             const ParticleProperties& def) const {
    if (has<ParametersList>(key)) {
      const auto& plist = get<ParametersList>(key);
      if (plist.has<pdgid_t>("pdgid") || plist.has<int>("pdgid"))
        return PDG::get()(plist.get<pdgid_t>("pdgid"));
      return ParticleProperties(plist);
    } else if (has<int>(key))
      return PDG::get()(get<int>(key));
    else {
      CG_DEBUG("ParametersList") << "Failed to retrieve particle properties parameter with key=" << key << ".";
      return def;
    }
  }

  /// Set a particle properties object value
  template <>
  ParametersList& ParametersList::set<ParticleProperties>(const std::string& key, const ParticleProperties& value) {
    return set<ParametersList>(key, value.parameters());
  }
}  // namespace cepgen

#undef IMPL_TYPE
