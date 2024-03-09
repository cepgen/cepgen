/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2018-2024  Laurent Forthomme
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

#include <iomanip>
#include <limits>
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
  }                                                                                                      \
  static_assert(true, "")

#define IMPL_TYPE_SET(type, coll, name)                                                                             \
  template <>                                                                                                       \
  bool ParametersList::has<type>(const std::string& key) const {                                                    \
    return coll.count(key) != 0;                                                                                    \
  }                                                                                                                 \
  template <>                                                                                                       \
  ParametersList& ParametersList::set<type>(std::string key, const type& value) {                                   \
    coll[std::move(key)] = (type)(value);                                                                           \
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
  }                                                                                                                 \
  template <>                                                                                                       \
  size_t ParametersList::erase<type>(const std::string& key) {                                                      \
    return coll.erase(key);                                                                                         \
  }                                                                                                                 \
  static_assert(true, "")

#define IMPL_TYPE_ALL(type, coll, name) \
  IMPL_TYPE_GET(type, coll, #name);     \
  IMPL_TYPE_SET(type, coll, #name);     \
  static_assert(true, "")

namespace cepgen {
  ParametersList::ParametersList(const ParametersList& oth) { operator+=(oth); }

  bool ParametersList::operator==(const ParametersList& oth) const {
    if (bool_values_ != oth.bool_values_)
      return false;
    if (int_values_ != oth.int_values_)
      return false;
    if (ulong_values_ != oth.ulong_values_)
      return false;
    if (dbl_values_ != oth.dbl_values_)
      return false;
    if (str_values_ != oth.str_values_)
      return false;
    if (lim_values_ != oth.lim_values_)
      return false;
    if (param_values_ != oth.param_values_)
      return false;
    if (vec_int_values_ != oth.vec_int_values_)
      return false;
    if (vec_dbl_values_ != oth.vec_dbl_values_)
      return false;
    if (vec_str_values_ != oth.vec_str_values_)
      return false;
    if (vec_lim_values_ != oth.vec_lim_values_)
      return false;
    if (vec_param_values_ != oth.vec_param_values_)
      return false;
    return true;
  }

  ParametersList& ParametersList::operator+=(const ParametersList& oth) {
    if (oth.empty() || *this == oth)  // ensure the two collections are not identical or empty
      return *this;
    if (empty()) {
      *this = oth;
      return *this;
    }
    // then check if any key of the other collection is already present in the list
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
      CG_DEBUG_LOOP("ParametersList") << utils::s("key", keys_erased.size(), true) << " erased: " << keys_erased << ".";
    //--- concatenate all typed lists
    bool_values_.insert(oth.bool_values_.begin(), oth.bool_values_.end());
    int_values_.insert(oth.int_values_.begin(), oth.int_values_.end());
    ulong_values_.insert(oth.ulong_values_.begin(), oth.ulong_values_.end());
    dbl_values_.insert(oth.dbl_values_.begin(), oth.dbl_values_.end());
    str_values_.insert(oth.str_values_.begin(), oth.str_values_.end());
    // special case for parameters collection: concatenate values instead of full containers
    for (const auto& par : oth.param_values_)
      // if the two parameters list are modules, and do not have the same name,
      // simply replace the old one with the new parameters list
      if (param_values_[par.first].getString(MODULE_NAME) == par.second.getString(MODULE_NAME))
        param_values_[par.first] += par.second;
      else
        param_values_[par.first] = par.second;
    lim_values_.insert(oth.lim_values_.begin(), oth.lim_values_.end());
    vec_int_values_.insert(oth.vec_int_values_.begin(), oth.vec_int_values_.end());
    vec_dbl_values_.insert(oth.vec_dbl_values_.begin(), oth.vec_dbl_values_.end());
    vec_str_values_.insert(oth.vec_str_values_.begin(), oth.vec_str_values_.end());
    vec_lim_values_.insert(oth.vec_lim_values_.begin(), oth.vec_lim_values_.end());
    vec_param_values_.insert(oth.vec_param_values_.begin(), oth.vec_param_values_.end());
    vec_vec_dbl_values_.insert(oth.vec_vec_dbl_values_.begin(), oth.vec_vec_dbl_values_.end());
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
    // first pre-process the arguments list to isolate all comma-separated arguments
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
                                    << " (split: " << utils::split(raw_list, ',') << "), "
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
      if (arg[arg.size() - 1] != '\'' && arg[arg.size() - 1] != '"' && cmd.size() > 1) {  // sub-parameters word found
        operator[]<ParametersList>(cmd.at(0)).feed(
            utils::merge(std::vector<std::string>(cmd.begin() + 1, cmd.end()), "/"));
        continue;
      }

      // from this moment on, a "key:value" or "key(:true)" was found
      const auto& subplist = utils::between(arg, "{", "}");
      if (!subplist.empty()) {
        for (const auto& subp : subplist)
          feed(subp);
        return *this;
      }
      const auto& word = cmd.at(0);
      auto words = utils::split(arg, ':');
      auto key = words.at(0);
      if (erase(key) > 0)
        CG_DEBUG("ParametersList:feed") << "Replacing key '" << key << "' with a new value.";
      if (key == "name")  // replace any "name" key encountered by the canonical module name key
        key = MODULE_NAME;
      if (words.size() == 1)  // basic key:true
        set<bool>(key, true);
      else if (words.size() == 2) {  // basic key:value
        const auto value = words.at(1);
        if (utils::isInt(value))
          set<int>(key, std::stoi(value));
        else if (utils::isFloat(value))
          set<double>(key, std::stod(value));
        else {
          const auto value_lc = utils::tolower(value);
          if (value_lc == "off" || value_lc == "no" || value_lc == "false")
            set<bool>(key, false);
          else if (value_lc == "on" || value_lc == "yes" || value_lc == "true")
            set<bool>(key, true);
          else if (value.find('>') != std::string::npos) {
            const auto limits = utils::split(value, '>');
            if (limits.size() != 2)
              throw CG_FATAL("ParametersList:feed") << "Failed to parse limits value '" << value << "'.";
            set<Limits>(key, Limits{std::stod(limits.at(0)), std::stod(limits.at(1))});
          } else {
            auto parsed_value = value;
            if (value.size() > 2 && value[0] == value[value.size() - 1] && (value[0] == '"' || value[0] == '\''))
              parsed_value = parsed_value.substr(1, value.size() - 2);
            set<std::string>(key, parsed_value);
          }
        }
      } else
        throw CG_FATAL("ParametersList:feed") << "Invalid key:value unpacking: " << word << "!";
    }
    return *this;
  }

  size_t ParametersList::erase(const std::string& key) {
    return erase<bool>(key) + erase<int>(key) + erase<unsigned long long>(key) + erase<double>(key) +
           erase<std::string>(key) + erase<Limits>(key) + erase<ParametersList>(key) + erase<std::vector<int> >(key) +
           erase<std::vector<double> >(key) + erase<std::vector<std::string> >(key) +
           erase<std::vector<ParametersList> >(key) + erase<std::vector<std::vector<double> > >(key);
  }

  bool ParametersList::empty() const { return keys(true).empty(); }

  std::ostream& operator<<(std::ostream& os, const ParametersList& params) {
    params.print(os);
    return os;
  }

  const ParametersList& ParametersList::print(std::ostream& os) const {
    const auto& keys_list = keys(true);
    if (keys_list.empty()) {
      os << "{}";
      return *this;
    }
    std::string sep;
    if (std::find(keys_list.begin(), keys_list.end(), MODULE_NAME) != keys_list.end()) {
      const auto plist_name = getNameString();
      auto mod_name = hasName<std::string>() ? "\"" + plist_name + "\"" : plist_name;
      os << "Module(" << mod_name, sep = ", ";
    } else
      os << "Parameters(";
    for (const auto& key : keys_list)
      if (key != MODULE_NAME)
        os << sep << key << "=" << getString(key, true), sep = ", ";
    os << ")";
    return *this;
  }

  std::string ParametersList::print() const {
    std::ostringstream os;
    print(os);
    return os.str();
  }

  std::vector<std::string> ParametersList::keys(bool name_key) const {
    std::vector<std::string> out{};
    auto key = [](const auto& p) { return p.first; };
    std::transform(bool_values_.begin(), bool_values_.end(), std::back_inserter(out), key);
    std::transform(int_values_.begin(), int_values_.end(), std::back_inserter(out), key);
    std::transform(ulong_values_.begin(), ulong_values_.end(), std::back_inserter(out), key);
    std::transform(dbl_values_.begin(), dbl_values_.end(), std::back_inserter(out), key);
    std::transform(str_values_.begin(), str_values_.end(), std::back_inserter(out), key);
    std::transform(param_values_.begin(), param_values_.end(), std::back_inserter(out), key);
    std::transform(lim_values_.begin(), lim_values_.end(), std::back_inserter(out), key);
    std::transform(vec_int_values_.begin(), vec_int_values_.end(), std::back_inserter(out), key);
    std::transform(vec_dbl_values_.begin(), vec_dbl_values_.end(), std::back_inserter(out), key);
    std::transform(vec_str_values_.begin(), vec_str_values_.end(), std::back_inserter(out), key);
    std::transform(vec_lim_values_.begin(), vec_lim_values_.end(), std::back_inserter(out), key);
    std::transform(vec_param_values_.begin(), vec_param_values_.end(), std::back_inserter(out), key);
    std::transform(vec_vec_dbl_values_.begin(), vec_vec_dbl_values_.end(), std::back_inserter(out), key);
    if (!name_key) {
      const auto it_name = std::find(out.begin(), out.end(), MODULE_NAME);
      if (it_name != out.end())
        out.erase(it_name);
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());  // at most 1 duplicate
    return out;
  }

  std::string ParametersList::getString(const std::string& key, bool wrap) const {
    // wrapper for the printout of a general variable
    auto wrap_val = [&wrap](const auto& val, const std::string& type) -> std::string {
      std::ostringstream os;
      if (type == "float" || type == "vfloat")
        os << std::fixed;
      else if (type == "bool")
        os << std::boolalpha;
      os << val;
      return (wrap ? type + "(" : "")  //+ (type == "bool" ? utils::yesno(std::stoi(os.str())) : os.str()) +
             + os.str() + (wrap ? ")" : "");
    };
    // wrapper for the printout of a collection type (vector, array, ...)
    auto wrap_coll = [&wrap_val](const auto& coll, const std::string& type) -> std::string {
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
    else if (has<std::vector<ParametersList> >(key))
      os << wrap_coll(get<std::vector<ParametersList> >(key), "VParams");
    else if (has<std::vector<int> >(key))
      os << wrap_coll(get<std::vector<int> >(key), "vint");
    else if (has<std::vector<double> >(key) || has<Limits>(key)) {
      std::string sep;
      if (has<std::vector<double> >(key)) {
        os << wrap_coll(get<std::vector<double> >(key), "vfloat");
        sep = "|";
      }
      if (has<Limits>(key))
        os << sep << wrap_val(get<Limits>(key), "Limits");
    } else if (has<std::vector<std::string> >(key))
      os << wrap_coll(get<std::vector<std::string> >(key), "vstr");
    else if (has<std::vector<Limits> >(key))
      os << wrap_coll(get<std::vector<Limits> >(key), "VLimits");
    else if (has<std::vector<std::vector<double> > >(key))
      os << wrap_coll(get<std::vector<std::vector<double> > >(key), "vvfloat");
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
    if (has<std::vector<Limits> >(old_key))
      set(new_key, get<std::vector<Limits> >(old_key)).erase(old_key);
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
        out << ":" << getString(key, false);
      sep = ",";
    }
    return out.str();
  }

  //------------------------------------------------------------------
  // default template (place holders)
  //------------------------------------------------------------------

  template <typename T>
  bool ParametersList::has(const std::string& key) const {
    throw CG_FATAL("ParametersList") << "Invalid type for key '" << key << "'.";
  }

  template <typename T>
  T ParametersList::get(const std::string& key, const T&) const {
    throw CG_FATAL("ParametersList") << "Invalid type retrieved for key '" << key << "'.";
  }

  template <typename T>
  T& ParametersList::operator[](const std::string& key) {
    throw CG_FATAL("ParametersList") << "Invalid type retrieved for key '" << key << "'.";
  }

  template <typename T>
  ParametersList& ParametersList::set(std::string key, const T&) {
    throw CG_FATAL("ParametersList") << "Invalid type to be set for key '" << key << "'.";
  }

  //------------------------------------------------------------------
  // sub-parameters-type attributes
  //------------------------------------------------------------------

  IMPL_TYPE_ALL(ParametersList, param_values_, "parameters");
  IMPL_TYPE_ALL(bool, bool_values_, "boolean");

  IMPL_TYPE_SET(int, int_values_, "integer");
  template <>
  int ParametersList::get<int>(const std::string& key, const int& def) const {
    if (has<int>(key))
      return int_values_.at(key);
    if (has<unsigned long long>(key)) {
      const auto ulong_val = ulong_values_.at(key);
      if (ulong_val >= std::numeric_limits<int>::max())
        CG_WARNING("ParametersList:get")
            << "Trying to retrieve a (too) long unsigned integer with an integer getter. Please fix your code.";
      return (int)ulong_val;
    }
    return def;
  }

  IMPL_TYPE_SET(unsigned long long, ulong_values_, "unsigned long integer");
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

  IMPL_TYPE_ALL(double, dbl_values_, "floating number");
  IMPL_TYPE_ALL(std::string, str_values_, "string");
  IMPL_TYPE_ALL(std::vector<int>, vec_int_values_, "vector of integers");
  IMPL_TYPE_ALL(std::vector<double>, vec_dbl_values_, "vector of floating numbers");
  IMPL_TYPE_ALL(std::vector<std::string>, vec_str_values_, "vector of strings");
  IMPL_TYPE_ALL(std::vector<Limits>, vec_lim_values_, "vector of limits");
  IMPL_TYPE_ALL(std::vector<ParametersList>, vec_param_values_, "vector of parameters");
  IMPL_TYPE_ALL(std::vector<std::vector<double> >, vec_vec_dbl_values_, "vector of vectors of floating numbers");

  //------------------------------------------------------------------
  // limits-type attributes
  //------------------------------------------------------------------

  IMPL_TYPE_SET(Limits, lim_values_, "limits");

  template <>
  Limits ParametersList::get<Limits>(const std::string& key, const Limits& def) const {
    // first try to find Limits object in collections
    Limits out;
    auto val = std::find_if(lim_values_.begin(), lim_values_.end(), [&key](const auto& kv) { return kv.first == key; });
    if (val != lim_values_.end())
      out = val->second;
    // still trying to build it from (min/max) attributes
    fill<double>(key + "min", out.min());
    fill<double>(key + "max", out.max());
    return out.validate();
    // nothing found ; returning default
    CG_DEBUG("ParametersList") << "Failed to retrieve limits parameter with key=" << key << ". "
                               << "Default value: " << def << ".";
    return def;
  }

  template <>
  const ParametersList& ParametersList::fill<Limits>(const std::string& key, Limits& value) const {
    fill<double>(key + "min", value.min());
    fill<double>(key + "max", value.max());
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

  //------------------------------------------------------------------
  // particle properties-type attributes
  //   particular case for this container, as it can either be
  //   represented by a ParametersList (collection of parameters) or
  //   an integer PDG identifier
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
    if (has<ParametersList>(key)) {  // first steer as a dictionary of particle properties
      const auto& plist = get<ParametersList>(key);
      if (plist.keys() == std::vector<std::string>{"pdgid"})
        return PDG::get()(plist.get<int>("pdgid"));
      return ParticleProperties(plist);
    } else if (has<pdgid_t>(key) ||
               has<int>(key)) {  // if not a dictionary of properties, retrieve from PDG runtime database
      CG_DEBUG("ParametersList") << "Retrieved physical properties for particle with PDG identifier '" << get<int>(key)
                                 << "' from PDG database.";
      return PDG::get()(get<int>(key));
    }
    CG_DEBUG("ParametersList") << "Failed to retrieve particle properties parameter with key=" << key << ".";
    return def;
  }

  /// Set a particle properties object value
  template <>
  ParametersList& ParametersList::set<ParticleProperties>(std::string key, const ParticleProperties& value) {
    PDG::get().define(value);
    return set<ParametersList>(std::move(key), value.parameters());
  }

  template <>
  std::vector<std::string> ParametersList::keysOf<ParticleProperties>() const {
    std::vector<std::string> pdesc_keys;
    for (const auto& key : keys())
      if (get<ParticleProperties>(key, ParticleProperties(-1)) != ParticleProperties(-1))
        pdesc_keys.emplace_back(key);
    return pdesc_keys;
  }
}  // namespace cepgen

#undef IMPL_TYPE_SET
#undef IMPL_TYPE_GET
#undef IMPL_TYPE_ALL
