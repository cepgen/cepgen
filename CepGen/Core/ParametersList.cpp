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

#include <iomanip>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/String.h"

#define IMPL_TYPE(type, coll, name)                                                                                 \
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
  type ParametersList::get<type>(const std::string& key, const type& def) const {                                   \
    auto val = std::find_if(coll.begin(), coll.end(), [&key](const auto& kv) { return kv.first == key; });          \
    if (val != coll.end())                                                                                          \
      return val->second;                                                                                           \
    CG_DEBUG("ParametersList") << "Failed to retrieve " << name << " parameter with key=" << key << ". "            \
                               << "Default value: " << def << ".";                                                  \
    return def;                                                                                                     \
  }                                                                                                                 \
  template <>                                                                                                       \
  std::vector<std::string> ParametersList::keysOf<type>() const {                                                   \
    std::vector<std::string> out;                                                                                   \
    std::transform(coll.begin(), coll.end(), std::back_inserter(out), [](const auto& pair) { return pair.first; }); \
    return out;                                                                                                     \
  }

namespace cepgen {
  const std::string ParametersList::MODULE_NAME = "mod_name";

  ParametersList::ParametersList(const ParametersList& oth)
      : param_values_(oth.param_values_),
        bool_values_(oth.bool_values_),
        int_values_(oth.int_values_),
        dbl_values_(oth.dbl_values_),
        str_values_(oth.str_values_),
        lim_values_(oth.lim_values_),
        vec_param_values_(oth.vec_param_values_),
        vec_int_values_(oth.vec_int_values_),
        vec_dbl_values_(oth.vec_dbl_values_),
        vec_str_values_(oth.vec_str_values_) {}

  bool ParametersList::operator==(const ParametersList& oth) const {
    if (keys() != oth.keys())
      return false;
    return true;
  }

  ParametersList& ParametersList::operator+=(const ParametersList& oth) {
    //--- first ensure no key is not already present in the list
    std::vector<std::string> keys_erased;
    for (const auto& key : oth.keys()) {
      if (has<ParametersList>(key)) {
        if (get<ParametersList>(key) == oth.get<ParametersList>(key) && erase(key) > 0)
          keys_erased.emplace_back(key);
      } else if (erase(key) > 0)
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
    lim_values_.insert(oth.lim_values_.begin(), oth.lim_values_.end());
    for (const auto& par : oth.param_values_)
      param_values_[par.first] += par.second;
    /*for (const auto& vpar : oth.vec_param_values_)
      for (const auto& par : vpar.second)
        if (!utils::contains(vec_param_values_[vpar.first], par))
          vec_param_values_[vpar.first].emplace_back(par);*/
    vec_param_values_.insert(oth.vec_param_values_.begin(), oth.vec_param_values_.end());
    return *this;
  }

  ParametersList ParametersList::operator+(const ParametersList& oth) const {
    ParametersList out = *this;
    out += oth;
    return out;
  }

  ParametersList& ParametersList::feed(const std::string& arg) {
    auto cmd = utils::split(arg, '/');
    auto& plist = cmd.size() > 1 ? operator[]<ParametersList>(cmd.at(0)) : *this;
    if (cmd.size() > 1)  // sub-parameters word found
      plist.feed(utils::merge(std::vector<std::string>(cmd.begin() + 1, cmd.end()), "/"));
    else {
      const auto word = cmd.at(0);
      auto words = utils::split(word, '=');
      auto key = words.at(0);
      if (key == "name")
        key = ParametersList::MODULE_NAME;
      if (words.size() == 1)  // basic key=true
        plist.set<bool>(key, true);
      else if (words.size() == 2) {  // basic key=value
        const auto value = words.at(1);
        try {
          if (value.find('.') != std::string::npos || value.find('e') != std::string::npos ||
              value.find('E') != std::string::npos)
            // double if contains a '.'/'e'
            plist.set<double>(key, std::stod(value));
          else
            plist.set<int>(key, std::stoi(value));
        } catch (const std::invalid_argument&) {
          const auto value_lc = utils::tolower(value);
          if (value_lc == "off" || value_lc == "no" || value_lc == "false")
            plist.set<bool>(key, false);
          else if (value_lc == "on" || value_lc == "yes" || value_lc == "true")
            plist.set<bool>(key, true);
          else
            plist.set<std::string>(key, value);
        }
      } else
        throw CG_FATAL("ParametersList:feed") << "Invalid key=value unpacking: " << word << "!";
    }
    return *this;
  }

  size_t ParametersList::erase(const std::string& key) {
    size_t out = 0ull;
    if (bool_values_.count(key) != 0)
      out += bool_values_.erase(key);
    if (int_values_.count(key) != 0)
      out += int_values_.erase(key);
    if (dbl_values_.count(key) != 0)
      out += dbl_values_.erase(key);
    if (str_values_.count(key) != 0)
      out += str_values_.erase(key);
    if (lim_values_.count(key) != 0)
      out += lim_values_.erase(key);
    if (param_values_.count(key) != 0)
      out += param_values_.erase(key);
    if (vec_int_values_.count(key) != 0)
      out += vec_int_values_.erase(key);
    if (vec_dbl_values_.count(key) != 0)
      out += vec_dbl_values_.erase(key);
    if (vec_str_values_.count(key) != 0)
      out += vec_str_values_.erase(key);
    if (vec_param_values_.count(key) != 0)
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
    const auto& mod_name = getString(ParametersList::MODULE_NAME);
    if (!mod_name.empty())
      os << "Module(" << mod_name, sep = ", ";
    else
      os << "Parameters(";
    for (const auto& key : keys(false))
      os << sep << key << "=" << getString(key, true), sep = ", ";
    os << ")";
    return *this;
  }

  std::vector<std::string> ParametersList::keys(bool name_key) const {
    std::vector<std::string> out;
    auto key = [](const auto& p) { return p.first; };
    std::transform(param_values_.begin(), param_values_.end(), std::back_inserter(out), key);
    std::transform(vec_param_values_.begin(), vec_param_values_.end(), std::back_inserter(out), key);
    std::transform(bool_values_.begin(), bool_values_.end(), std::back_inserter(out), key);
    std::transform(int_values_.begin(), int_values_.end(), std::back_inserter(out), key);
    std::transform(vec_int_values_.begin(), vec_int_values_.end(), std::back_inserter(out), key);
    std::transform(dbl_values_.begin(), dbl_values_.end(), std::back_inserter(out), key);
    std::transform(vec_dbl_values_.begin(), vec_dbl_values_.end(), std::back_inserter(out), key);
    std::transform(str_values_.begin(), str_values_.end(), std::back_inserter(out), key);
    std::transform(vec_str_values_.begin(), vec_str_values_.end(), std::back_inserter(out), key);
    std::transform(lim_values_.begin(), lim_values_.end(), std::back_inserter(out), key);
    if (!name_key) {
      const auto it_name = std::find(out.begin(), out.end(), MODULE_NAME);
      if (it_name != out.end())
        out.erase(it_name);
    }
    std::sort(out.begin(), out.end());
    return out;
  }

  std::string ParametersList::getString(const std::string& key, bool wrap) const {
    std::ostringstream os;
    auto wrap_val = [&wrap](const auto& val, const std::string& type) -> std::string {
      std::ostringstream os;
      os << val;
      return (wrap ? type + "(" : "") + (type == "bool" ? utils::yesno(std::stoi(os.str())) : os.str()) +
             (wrap ? ")" : "");
    };
    auto wrap_coll = [&wrap, &wrap_val](const auto& coll, const std::string& type) -> std::string {
      return wrap_val(utils::merge(coll, ", "), type);
    };
    if (has<ParametersList>(key)) {
      auto plist = get<ParametersList>(key);
      if (!wrap)
        os << plist;
      else {
        const auto& plist_name = plist.getString(MODULE_NAME, false);
        if (plist_name.empty())
          os << "Parameters(" << plist << ")";
        else {
          plist.erase(MODULE_NAME);
          os << "Module(" << plist_name << ", " << plist << ")";
        }
      }
    } else if (has<bool>(key))
      os << wrap_val(get<bool>(key), "bool");
    else if (has<int>(key))
      os << wrap_val(get<int>(key), "int");
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

  IMPL_TYPE(bool, bool_values_, "boolean")
  IMPL_TYPE(int, int_values_, "integer")
  IMPL_TYPE(std::vector<int>, vec_int_values_, "vector of integers")
  IMPL_TYPE(double, dbl_values_, "floating number")
  IMPL_TYPE(std::vector<double>, vec_dbl_values_, "vector of floating numbers")
  IMPL_TYPE(std::string, str_values_, "string")
  IMPL_TYPE(std::vector<std::string>, vec_str_values_, "vector of strings")
  IMPL_TYPE(ParametersList, param_values_, "parameters")
  IMPL_TYPE(std::vector<ParametersList>, vec_param_values_, "vector of parameters")

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
    return param_values_.count(key) != 0;
  }

  /// Get a particle properties object
  template <>
  ParticleProperties ParametersList::get<ParticleProperties>(const std::string& key,
                                                             const ParticleProperties& def) const {
    if (has<ParametersList>(key)) {
      const auto& plist = get<ParametersList>(key);
      ParticleProperties out;
      const auto& pdgid = plist.get<int>("pdgid", 0);
      if (PDG::get().has(pdgid))
        out = PDG::get()(pdgid);
      else
        out.pdgid = pdgid;
      bool modified = false;
      if (plist.has<std::string>("name"))
        out.name = plist.get<std::string>("name"), modified = true;
      if (plist.has<std::string>("description"))
        out.description = plist.get<std::string>("description"), modified = true;
      if (plist.has<int>("colours"))
        out.colours = plist.get<int>("colours"), modified = true;
      if (plist.has<double>("mass"))
        out.mass = plist.get<double>("mass"), modified = true;
      if (plist.has<double>("width"))
        out.width = plist.get<double>("width"), modified = true;
      if (plist.has<double>("charge"))
        out.charge = short(plist.get<double>("charge") * 3), modified = true;
      if (plist.has<bool>("fermion"))
        out.fermion = plist.get<bool>("fermion"), modified = true;
      if (modified)
        PDG::get().define(out);
      return out;
    } else if (has<int>(key))
      return PDG::get()(get<int>(key));
    else {
      CG_DEBUG("ParametersList") << "Failed to retrieve parameter with key=" << key << ".";
      return def;
    }
  }

  /// Set a particle properties object value
  template <>
  ParametersList& ParametersList::set<ParticleProperties>(const std::string& key, const ParticleProperties& value) {
    return set<ParametersList>(key,
                               ParametersList()
                                   .set<int>("pdgid", value.pdgid)
                                   .set<std::string>("name", value.name)
                                   .set<std::string>("description", value.description)
                                   .set<int>("colours", value.colours)
                                   .set<double>("mass", value.mass)
                                   .set<double>("width", value.width)
                                   .set<double>("charge", value.charge * 1. / 3)
                                   .set<bool>("fermion", value.fermion));
  }
}  // namespace cepgen

#undef IMPL_TYPE
