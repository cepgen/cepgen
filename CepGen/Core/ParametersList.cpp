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

namespace cepgen {
  const std::string ParametersList::MODULE_NAME = "mod_name";

  ParametersList::ParametersList(const ParametersList& oth)
      : param_values_(oth.param_values_),
        int_values_(oth.int_values_),
        dbl_values_(oth.dbl_values_),
        str_values_(oth.str_values_),
        lim_values_(oth.lim_values_),
        vec_param_values_(oth.vec_param_values_),
        vec_int_values_(oth.vec_int_values_),
        vec_dbl_values_(oth.vec_dbl_values_),
        vec_str_values_(oth.vec_str_values_) {}

  ParametersList& ParametersList::operator+=(const ParametersList& oth) {
    //--- first ensure no key is not already present in the list
    for (const auto& key : oth.keys())
      erase(key);
    //--- concatenate all typed lists
    int_values_.insert(oth.int_values_.begin(), oth.int_values_.end());
    dbl_values_.insert(oth.dbl_values_.begin(), oth.dbl_values_.end());
    str_values_.insert(oth.str_values_.begin(), oth.str_values_.end());
    lim_values_.insert(oth.lim_values_.begin(), oth.lim_values_.end());
    vec_param_values_.insert(oth.vec_param_values_.begin(), oth.vec_param_values_.end());
    vec_int_values_.insert(oth.vec_int_values_.begin(), oth.vec_int_values_.end());
    vec_dbl_values_.insert(oth.vec_dbl_values_.begin(), oth.vec_dbl_values_.end());
    vec_str_values_.insert(oth.vec_str_values_.begin(), oth.vec_str_values_.end());
    for (const auto& par : oth.param_values_)
      param_values_[par.first] += par.second;
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
            plist.set<int>(key, std::stod(value));
        } catch (const std::invalid_argument&) {
          if (value == "off" || value == "no" || value == "false")
            plist.set<bool>(key, false);
          else if (value == "on" || value == "yes" || value == "true")
            plist.set<bool>(key, true);
          else
            plist.set<std::string>(key, value);
        }
      } else
        throw CG_FATAL("ParametersList:feed") << "Invalid key=value unpacking: " << word << "!";
    }
    return *this;
  }

  size_t ParametersList::erase(std::string key) {
    size_t out = 0ull;
    if (int_values_.count(key) != 0)
      out += int_values_.erase(key);
    if (dbl_values_.count(key) != 0)
      out += dbl_values_.erase(key);
    if (str_values_.count(key) != 0)
      out += str_values_.erase(key);
    if (lim_values_.count(key) != 0)
      out += lim_values_.erase(key);
    if (vec_param_values_.count(key) != 0)
      out += vec_param_values_.erase(key);
    if (vec_int_values_.count(key) != 0)
      out += vec_int_values_.erase(key);
    if (vec_dbl_values_.count(key) != 0)
      out += vec_dbl_values_.erase(key);
    if (vec_str_values_.count(key) != 0)
      out += vec_str_values_.erase(key);
    if (param_values_.count(key) != 0)
      out += param_values_.erase(key);
    return out;
  }

  bool ParametersList::empty() const { return keys(false).empty(); }

  std::ostream& operator<<(std::ostream& os, const ParametersList& params) {
    std::string sep;
    const auto key_name = [&](const auto& key) -> std::string {
      return key == ParametersList::MODULE_NAME ? "name=" : key + "=";
    };
    const auto mod_or_param = [&](const auto& plist) {
      const auto& pkeys = plist.keys();
      auto name = std::find(pkeys.begin(), pkeys.end(), ParametersList::MODULE_NAME);
      std::ostringstream os;
      if (name != pkeys.end()) {
        auto plist_copy = plist;
        plist_copy.erase(ParametersList::MODULE_NAME);
        os << "Module(" << plist.getString(ParametersList::MODULE_NAME) << ", " << plist_copy << ")";
      } else
        os << "Parameters(" << plist << ")";
      return os.str();
    };
    for (const auto& kv : params.int_values_)
      os << sep << key_name(kv.first) << "int(" << kv.second << ")", sep = ", ";
    for (const auto& kv : params.dbl_values_)
      os << sep << key_name(kv.first) << "float(" << kv.second << ")", sep = ", ";
    for (const auto& kv : params.str_values_)
      os << sep << key_name(kv.first) << "str(" << kv.second << ")", sep = ", ";
    for (const auto& kv : params.param_values_)
      os << sep << key_name(kv.first) << mod_or_param(kv.second), sep = ", ";
    for (const auto& kv : params.lim_values_)
      os << sep << key_name(kv.first) << "limits(" << kv.second << ")", sep = ", ";
    for (const auto& kv : params.vec_int_values_) {
      os << sep << key_name(kv.first) << "vint(", sep = ", ";
      std::string sep1;
      for (const auto& val : kv.second)
        os << sep1 << val, sep1 = ", ";
      os << ")";
    }
    for (const auto& kv : params.vec_dbl_values_) {
      os << sep << key_name(kv.first) << "vfloat(", sep = ", ";
      std::string sep1;
      for (const auto& val : kv.second)
        os << sep1 << val, sep1 = ", ";
      os << ")";
    }
    for (const auto& kv : params.vec_str_values_) {
      os << sep << key_name(kv.first) << "vstr(", sep = ", ";
      std::string sep1;
      for (const auto& val : kv.second)
        os << sep1 << val, sep1 = ", ";
      os << ")";
    }
    return os;
  }

  std::vector<std::string> ParametersList::keys(bool name_key) const {
    std::vector<std::string> out;
    auto key = [](const auto& p) { return p.first; };
    std::transform(param_values_.begin(), param_values_.end(), std::back_inserter(out), key);
    std::transform(vec_param_values_.begin(), vec_param_values_.end(), std::back_inserter(out), key);
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
    return out;
  }

  std::string ParametersList::getString(const std::string& key) const {
    std::ostringstream os;
    if (has<ParametersList>(key))
      os << "params{" << get<ParametersList>(key) << "}";
    else if (has<int>(key))
      os << get<int>(key);
    else if (has<double>(key))
      os << get<double>(key);
    else if (has<std::string>(key))
      os << get<std::string>(key);
    else if (has<Limits>(key))
      os << get<Limits>(key);
    else if (has<std::vector<ParametersList> >(key)) {
      std::string sep;
      for (const auto& p : get<std::vector<ParametersList> >(key))
        os << sep << p, sep = ", ";
    } else if (has<std::vector<int> >(key)) {
      std::string sep;
      for (const auto& p : get<std::vector<int> >(key))
        os << sep << p, sep = ", ";
    } else if (has<std::vector<double> >(key)) {
      std::string sep;
      for (const auto& p : get<std::vector<double> >(key))
        os << sep << p, sep = ", ";
    } else if (has<std::vector<std::string> >(key)) {
      std::string sep;
      for (const auto& p : get<std::vector<std::string> >(key))
        os << sep << p, sep = ", ";
    }
    return os.str();
  }

  //------------------------------------------------------------------
  // default template (placeholders)
  //------------------------------------------------------------------

  template <typename T>
  bool ParametersList::has(std::string key) const {
    throw CG_FATAL("ParametersList") << "Invalid type for key=" << key << "!";
  }

  template <typename T>
  T ParametersList::get(std::string key, const T& def) const {
    throw CG_FATAL("ParametersList") << "Invalid type retrieved for key=" << key << "!";
  }

  template <typename T>
  T& ParametersList::operator[](std::string key) {
    throw CG_FATAL("ParametersList") << "Invalid type retrieved for key=" << key << "!";
  }

  template <typename T>
  ParametersList& ParametersList::set(std::string key, const T& value) {
    throw CG_FATAL("ParametersList") << "Invalid type to be set for key=" << key << "!";
  }

  //------------------------------------------------------------------
  // sub-parameters-type attributes
  //------------------------------------------------------------------

  template <>
  ParametersList ParametersList::get<ParametersList>(std::string key, const ParametersList& def) const {
    auto val =
        std::find_if(param_values_.begin(), param_values_.end(), [&key](const auto& kv) { return kv.first == key; });
    if (val != param_values_.end())
      return val->second;
    CG_DEBUG("ParametersList") << "Failed to retrieve parameters with key=" << key << ". "
                               << "Default value: " << def << ".";
    return def;
  }

  template <>
  std::vector<ParametersList> ParametersList::get<std::vector<ParametersList> >(
      std::string key, const std::vector<ParametersList>& def) const {
    auto val = std::find_if(
        vec_param_values_.begin(), vec_param_values_.end(), [&key](const auto& kv) { return kv.first == key; });
    if (val != vec_param_values_.end())
      return val->second;
    CG_DEBUG("ParametersList") << "Failed to retrieve parameters collection with key=" << key << ". "
                               << "Default value: " << def << ".";
    return def;
  }

  //------------------------------------------------------------------
  // integer-type attributes
  //------------------------------------------------------------------

  template <>
  int ParametersList::get<int>(std::string key, const int& def) const {
    auto val = std::find_if(int_values_.begin(), int_values_.end(), [&key](const auto& kv) { return kv.first == key; });
    if (val != int_values_.end())
      return val->second;
    CG_DEBUG("ParametersList") << "Failed to retrieve integer parameter with key=" << key << ". "
                               << "Default value: " << def << ".";
    return def;
  }

  template <>
  std::vector<int> ParametersList::get<std::vector<int> >(std::string key, const std::vector<int>& def) const {
    auto val = std::find_if(
        vec_int_values_.begin(), vec_int_values_.end(), [&key](const auto& kv) { return kv.first == key; });
    if (val != vec_int_values_.end())
      return val->second;
    CG_DEBUG("ParametersList") << "Failed to retrieve integer collection with key=" << key << ". "
                               << "Default value: " << def << ".";
    return def;
  }

  //------------------------------------------------------------------
  // floating point-type attributes
  //------------------------------------------------------------------

  template <>
  double ParametersList::get<double>(std::string key, const double& def) const {
    auto val = std::find_if(dbl_values_.begin(), dbl_values_.end(), [&key](const auto& kv) { return kv.first == key; });
    if (val != dbl_values_.end())
      return val->second;
    CG_DEBUG("ParametersList") << "Failed to retrieve double parameter with key=" << key << ". "
                               << "Default value: " << def << ".";
    return def;
  }

  template <>
  std::vector<double> ParametersList::get<std::vector<double> >(std::string key, const std::vector<double>& def) const {
    auto val = std::find_if(
        vec_dbl_values_.begin(), vec_dbl_values_.end(), [&key](const auto& kv) { return kv.first == key; });
    if (val != vec_dbl_values_.end())
      return val->second;
    CG_DEBUG("ParametersList") << "Failed to retrieve double collection with key=" << key << ". "
                               << "Default value: " << def << ".";
    return def;
  }

  //------------------------------------------------------------------
  // string-type attributes
  //------------------------------------------------------------------

  template <>
  std::string ParametersList::get<std::string>(std::string key, const std::string& def) const {
    auto val = std::find_if(str_values_.begin(), str_values_.end(), [&key](const auto& kv) { return kv.first == key; });
    if (val != str_values_.end())
      return val->second;
    CG_DEBUG("ParametersList") << "Failed to retrieve string parameter with key=" << key << ". "
                               << "Default value: " << def << ".";
    return def;
  }

  template <>
  std::vector<std::string> ParametersList::get<std::vector<std::string> >(std::string key,
                                                                          const std::vector<std::string>& def) const {
    auto val = std::find_if(
        vec_str_values_.begin(), vec_str_values_.end(), [&key](const auto& kv) { return kv.first == key; });
    if (val != vec_str_values_.end())
      return val->second;
    CG_DEBUG("ParametersList") << "Failed to retrieve string collection with key=" << key << ". "
                               << "Default value: " << def << ".";
    return def;
  }

  //------------------------------------------------------------------
  // limits-type attributes
  //------------------------------------------------------------------

  template <>
  Limits ParametersList::get<Limits>(std::string key, const Limits& def) const {
    // first try to find Limits object in collections
    auto val = std::find_if(lim_values_.begin(), lim_values_.end(), [&key](const auto& kv) { return kv.first == key; });
    if (val != lim_values_.end())
      return val->second;
    // Limits object not found ; still trying to build it from (min/max) attributes
    Limits buf;
    if (has<double>(key + "min"))
      fill<double>(key + "min", buf.min());
    if (has<double>(key + "max"))
      fill<double>(key + "max", buf.max());
    if (buf.valid())
      return buf.validate();
    // nothing found ; returning default
    CG_DEBUG("ParametersList") << "Failed to retrieve limits parameter with key=" << key << ". "
                               << "Default value: " << def << ".";
    return def;
  }

  template <>
  const ParametersList& ParametersList::fill<Limits>(std::string key, Limits& value) const {
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

  //------------------------------------------------------------------
  // particle properties-type attributes
  //------------------------------------------------------------------

  template <>
  ParticleProperties ParametersList::get<ParticleProperties>(std::string key, const ParticleProperties& def) const {
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

  template <>
  ParametersList& ParametersList::set<ParticleProperties>(std::string key, const ParticleProperties& value) {
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
