/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2018-2025  Laurent Forthomme
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
#include <regex>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/String.h"

using namespace cepgen;
using namespace std::string_literals;

#define IMPL_TYPE_GET(type, values_collection)                                                                  \
  template <>                                                                                                   \
  type ParametersList::get<type>(const std::string& key, const type& default_value) const {                     \
    if (values_collection.count(key) > 0)                                                                       \
      return values_collection.at(key);                                                                         \
    CG_DEBUG("ParametersList") << "Failed to retrieve " << utils::demangle(typeid(type).name())                 \
                               << " parameter with key=" << key << ". Default value: " << default_value << "."; \
    return default_value;                                                                                       \
  }                                                                                                             \
  static_assert(true, "")

#define IMPL_TYPE_SET(type, values_collection)                                                                         \
  template <>                                                                                                          \
  bool ParametersList::has<type>(const std::string& key) const {                                                       \
    return values_collection.count(key) != 0;                                                                          \
  }                                                                                                                    \
  template <>                                                                                                          \
  ParametersList& ParametersList::set<type>(const std::string& key, const type& value) {                               \
    values_collection[key] = static_cast<type>(value);                                                                 \
    return *this;                                                                                                      \
  }                                                                                                                    \
  template <>                                                                                                          \
  type& ParametersList::operator[]<type>(const std::string& key) {                                                     \
    return values_collection[key];                                                                                     \
  }                                                                                                                    \
  template <>                                                                                                          \
  std::vector<std::string> ParametersList::keysOf<type>() const {                                                      \
    std::vector<std::string> out;                                                                                      \
    std::transform(values_collection.begin(), values_collection.end(), std::back_inserter(out), [](const auto& pair) { \
      return pair.first;                                                                                               \
    });                                                                                                                \
    return out;                                                                                                        \
  }                                                                                                                    \
  template <>                                                                                                          \
  size_t ParametersList::erase<type>(const std::string& key) {                                                         \
    return values_collection.erase(key);                                                                               \
  }                                                                                                                    \
  static_assert(true, "")

#define IMPL_TYPE_ALL(type, values_collection) \
  IMPL_TYPE_GET(type, values_collection);      \
  IMPL_TYPE_SET(type, values_collection);      \
  static_assert(true, "")

ParametersList::ParametersList(const ParametersList& oth) { operator+=(oth); }

bool ParametersList::operator==(const ParametersList& oth) const {
#define __TYPE_ENUM(type, map, name) \
  if (map != oth.map)                \
    return false;
  REGISTER_CONTENT_TYPE
#undef __TYPE_ENUM
  return true;
}

ParametersList ParametersList::diff(const ParametersList& oth) const {
  ParametersList diff;
  ParametersList &mine = diff.operator[]<ParametersList>("mine"), &theirs = diff.operator[]<ParametersList>("theirs");
  if (*this == oth)
    return diff;
  for (const auto& key : keys()) {
    if (has<ParametersList>(key)) {
      const auto& my_plist = get<ParametersList>(key);
      if (!has<ParametersList>(key))
        mine.set(key, my_plist);
      else if (const auto their_plist = oth.get<ParametersList>(key); my_plist != their_plist) {
        mine.set(key, my_plist);
        theirs.set(key, their_plist);
      }
      continue;
    }
#define __TYPE_ENUM(type, map, name)                                                                     \
  if (const auto my_param = get<type>(key), their_param = oth.get<type>(key); my_param != their_param) { \
    mine.set(key, my_param);                                                                             \
    if (!oth.empty())                                                                                    \
      theirs.set(key, their_param);                                                                      \
    continue;                                                                                            \
  }
    REGISTER_CONTENT_TYPE
#undef __TYPE_ENUM
  }
  return diff;
}

ParametersList& ParametersList::operator+=(const ParametersList& other) {
  if (empty()) {
    *this = other;
    return *this;
  }
  if (other.empty() || *this == other)  // ensure the two collections are not identical or empty
    return *this;
  auto other_modified = other;
  std::vector<std::string> keys_erased;
  for (const auto& other_key :
       other_modified.keys()) {  // check if any key of the other collection is already present in the list
    if (has<ParametersList>(other_key)) {
      // do not remove a duplicate parameters collection if they are not strictly identical ;
      // will concatenate its values with the other object's
      if (get<ParametersList>(other_key) == other_modified.get<ParametersList>(other_key) &&
          other_modified.erase(other_key) > 0)
        keys_erased.emplace_back(other_key);
    } else if (other_modified.has<double>(other_key) &&
               (utils::endsWith(other_key, "max") || utils::endsWith(other_key, "min"))) {
      // hacky path to drop '*min'/'*max' keys if '*' limits were found
      if (auto& limits = operator[]<Limits>(other_key.substr(0, other_key.size() - 3));
          utils::endsWith(other_key, "max"))
        limits.max() = other_modified.get<double>(other_key);
      else if (utils::endsWith(other_key, "min"))
        limits.min() = other_modified.get<double>(other_key);
      if (erase(other_key) > 0 && other_modified.erase(other_key) > 0)
        keys_erased.emplace_back(other_key);
    } else if (erase(other_key) > 0)  // any other duplicate key is just replaced
      keys_erased.emplace_back(other_key);
  }
  if (!keys_erased.empty())
    CG_DEBUG_LOOP("ParametersList") << utils::s("key", keys_erased.size(), true) << " erased: " << keys_erased << ".";
  // concatenate all typed lists
#define __TYPE_ENUM(type, map, name) map.insert(other_modified.map.begin(), other_modified.map.end());
  REGISTER_CONTENT_TYPE
#undef __TYPE_ENUM
  // special case for parameters collection: concatenate values instead of full containers
  for (const auto& [key, parameters] : other_modified.param_values_)
    // if the two parameters list are modules, and do not have the same name,
    // simply replace the old one with the new parameters list
    if (param_values_[key].getNameString() == parameters.getNameString())
      param_values_[key] += parameters;
    else
      param_values_[key] = parameters;
  return *this;
}

ParametersList ParametersList::operator+(const ParametersList& other) const {
  auto out = *this;
  out += other;
  return out;
}

ParametersList& ParametersList::feed(const std::string& raw_args) {
  auto raw_list = raw_args;
  if (const auto raw_list_stripped = utils::between(raw_list, "{", "}");
      raw_list_stripped.size() == 1 && raw_list == "{" + raw_list_stripped[0] + "}")
    raw_list = raw_list_stripped[0];
  // first pre-process the arguments list to isolate all comma-separated arguments
  std::vector<std::string> arguments_list;
  std::vector<std::string> buffer;
  short num_open_braces = 0;
  for (const auto& item : utils::split(raw_list, ',')) {
    buffer.emplace_back(item);
    num_open_braces += std::count(item.begin(), item.end(), '{') - std::count(item.begin(), item.end(), '}');
    if (num_open_braces <= 0) {
      arguments_list.emplace_back(utils::merge(buffer, ","));
      buffer.clear();
    }
  }
  CG_DEBUG("ParametersList:feed") << "Parsed arguments: " << arguments_list << ", raw list: " << raw_list
                                  << " (split: " << utils::split(raw_list, ',') << "), "
                                  << "{-} imbalance: " << num_open_braces << ".";
  if (num_open_braces != 0)
    throw CG_ERROR("ParametersList:feed") << "Invalid string to be parsed as a parameters list!\n\t"
                                          << "Open-closed braces imbalance: " << num_open_braces << "\n\t"
                                          << "Raw list: " << raw_list << "\n\t"
                                          << "Resulting list: " << arguments_list << ", buffer: " << buffer << ".";
  for (const auto& argument : arguments_list) {  // loop through all unpacked arguments
    auto command = utils::split(argument, '/');  // browse through the parameters hierarchy
    if (argument[argument.size() - 1] != '\'' && argument[argument.size() - 1] != '"' &&
        command.size() > 1) {  // sub-parameters word found
      operator[]<ParametersList>(command.at(0)).feed(utils::merge(std::vector(command.begin() + 1, command.end()), "/"));
      continue;
    }
    // from this point, a "key:value" or "key(:true)" was found
    if (const auto& sub_plist = utils::between(argument, "{", "}"); !sub_plist.empty()) {
      for (const auto& subp : sub_plist)
        feed(subp);
      return *this;
    }
    const auto& word = command.at(0);
    const auto words = utils::split(argument, ':');
    auto key = words.at(0);
    if (erase(key) > 0)
      CG_DEBUG("ParametersList:feed") << "Replacing key '" << key << "' with a new value.";
    if (key == "name")  // replace any "name" key encountered by the canonical module name key
      key = MODULE_NAME;
    if (words.size() == 1)  // basic key:true
      set<bool>(key, true);
    else if (words.size() == 2) {  // basic key:value
      if (const auto& value = words.at(1); utils::isInt(value))
        set(key, std::stoi(value));
      else if (utils::isFloat(value))
        set(key, std::stod(value));
      else if (value[0] == value[value.size() - 1] && (value[0] == '\'' || value[0] == '"'))  // string parameter
        set(key, value.substr(1, value.size() - 2));
      else {
        if (const auto value_lc = utils::toLower(value); value_lc == "off" || value_lc == "no" || value_lc == "false")
          set<bool>(key, false);
        else if (value_lc == "on" || value_lc == "yes" || value_lc == "true")
          set<bool>(key, true);
        else if (value.find('>') != std::string::npos) {
          const auto limits = utils::split(value, '>');
          if (limits.size() != 2)
            throw CG_FATAL("ParametersList:feed") << "Failed to parse limits value '" << value << "'.";
          set(key, Limits{std::stod(limits.at(0)), std::stod(limits.at(1))});
        } else {
          auto parsed_value = value;
          if (value.size() > 2 && value[0] == value[value.size() - 1] && (value[0] == '"' || value[0] == '\''))
            parsed_value = parsed_value.substr(1, value.size() - 2);
          set(key, parsed_value);
        }
      }
    } else
      throw CG_FATAL("ParametersList:feed") << "Invalid key:value unpacking: " << word << "!";
  }
  return *this;
}

size_t ParametersList::erase(const std::string& key) {
  size_t num_keys_erased = 0;
#define __TYPE_ENUM(type, map, name) num_keys_erased += erase<type>(key);
  REGISTER_CONTENT_TYPE
#undef __TYPE_ENUM
  return num_keys_erased;
}

bool ParametersList::empty() const { return keys(true).empty(); }

const ParametersList& ParametersList::print(std::ostream& os) const {
  const auto& keys_list = keys(true);
  if (keys_list.empty()) {
    os << "{}";
    return *this;
  }
  std::string sep;
  if (std::find(keys_list.begin(), keys_list.end(), MODULE_NAME) != keys_list.end()) {
    const auto plist_name = getNameString();
    const auto mod_name = hasName() ? "\"" + plist_name + "\"" : plist_name;
    os << "Module(" << mod_name, sep = ", ";
  } else
    os << "Parameters(";
  for (const auto& key : keys_list)
    if (key != MODULE_NAME)
      os << sep << key << "=" << getString(key, true), sep = ", ";
  os << ")";
  return *this;
}

std::string ParametersList::print(bool compact) const {
  std::ostringstream os;
  if (compact) {
    os << utils::colourise(getNameString(), utils::Colour::none, utils::Modifier::underline);
    if (const auto& keys_list = keys(false); !keys_list.empty()) {
      std::string sep = "{";
      for (const auto& key : keys_list)
        os << sep << key << "=" << getString(key), sep = ", ";
      os << "}";
    }
    return os.str();
  }
  print(os);
  return os.str();
}

bool ParametersList::hasName() const { return has<std::string>(MODULE_NAME); }

std::string ParametersList::name(const std::string& default_value) const {
  return get<std::string>(MODULE_NAME, default_value);
}

ParametersList& ParametersList::setName(const std::string& value) { return set(MODULE_NAME, value); }

std::vector<std::string> ParametersList::keys(bool name_key) const {
  std::vector<std::string> out{};
  const auto key = [](const auto& p) { return p.first; };
#define __TYPE_ENUM(type, map, name) std::transform(map.begin(), map.end(), std::back_inserter(out), key);
  REGISTER_CONTENT_TYPE
#undef __TYPE_ENUM
  if (!name_key) {
    if (const auto it_name = std::find(out.begin(), out.end(), MODULE_NAME); it_name != out.end())
      out.erase(it_name);
  }
  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());  // at most 1 duplicate
  return out;
}

std::string ParametersList::getString(const std::string& key, bool wrap) const {
  // wrapper for the printout of a general variable
  auto wrap_val = [&wrap](const auto& val, const std::string& type) -> std::string {
    return (wrap ? type + "(" : "") + utils::merge(val, ",") + (wrap ? ")" : "");
  };
  // wrapper for the printout of a collection type (vector, array, ...)
  auto wrap_coll = [&wrap_val](const auto& coll, const std::string& type) -> std::string {
    return wrap_val(utils::merge(coll, ", "), type);
  };
  std::ostringstream os;
  if (has<ParametersList>(key)) {
    os << get<ParametersList>(key);
    return os.str();
  }
  if (has<std::vector<double> >(key) || has<Limits>(key)) {
    std::string sep;
    if (has<std::vector<double> >(key)) {
      os << wrap_coll(get<std::vector<double> >(key), "vfloat"s);
      sep = "|";
    }
    if (has<Limits>(key))
      os << sep << wrap_val(get<Limits>(key), "Limits");
    return os.str();
  }
  if (has<bool>(key)) {
    os << std::boolalpha << get<bool>(key);
    return os.str();
  }
#define __TYPE_ENUM(type, map, name) \
  if (has<type>(key))                \
    return wrap_val(get<type>(key), name);
  REGISTER_CONTENT_TYPE
#undef __TYPE_ENUM
  if (key == MODULE_NAME)
    return "";
  throw CG_ERROR("ParametersList:getString")
      << "Unrecognised type for key '" << key << "' from parameters list " << *this << ".";
}

ParametersList& ParametersList::rename(const std::string& old_key, const std::string& new_key) {
#define __TYPE_ENUM(type, map, name) \
  if (has<type>(old_key))            \
    set(new_key, get<type>(old_key)).erase(old_key);
  REGISTER_CONTENT_TYPE
#undef __TYPE_ENUM
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
size_t ParametersList::erase(const std::string& key) {
  throw CG_FATAL("ParametersList") << "Invalid type for key '" << key << "'.";
}

template <typename T>
bool ParametersList::has(const std::string& key) const {
  throw CG_FATAL("ParametersList") << "Invalid type for key '" << key << "'.";
}

template <typename T>
std::vector<std::string> ParametersList::keysOf() const {
  throw CG_FATAL("ParametersList") << "Invalid type to retrieve keys from.";
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
ParametersList& ParametersList::set(const std::string& key, const T&) {
  throw CG_FATAL("ParametersList") << "Invalid type to be set for key '" << key << "'.";
}

//------------------------------------------------------------------
// sub-parameters-type attributes
//------------------------------------------------------------------

IMPL_TYPE_ALL(ParametersList, param_values_);
IMPL_TYPE_ALL(bool, bool_values_);
IMPL_TYPE_ALL(int, int_values_);
IMPL_TYPE_ALL(unsigned long long, ulong_values_);
IMPL_TYPE_ALL(double, dbl_values_);
IMPL_TYPE_ALL(std::string, str_values_);
IMPL_TYPE_ALL(std::vector<int>, vec_int_values_);
IMPL_TYPE_ALL(std::vector<double>, vec_dbl_values_);
IMPL_TYPE_ALL(std::vector<std::string>, vec_str_values_);
IMPL_TYPE_ALL(std::vector<Limits>, vec_lim_values_);
IMPL_TYPE_ALL(std::vector<ParametersList>, vec_param_values_);
IMPL_TYPE_ALL(std::vector<std::vector<double> >, vec_vec_dbl_values_);

//------------------------------------------------------------------
// limits-type attributes
//------------------------------------------------------------------

IMPL_TYPE_SET(Limits, lim_values_);

template <>
Limits ParametersList::get<Limits>(const std::string& key, const Limits& default_value) const {
  // first try to find Limits object in collections
  auto out = default_value;
  if (const auto val =
          std::find_if(lim_values_.begin(), lim_values_.end(), [&key](const auto& kv) { return kv.first == key; });
      val != lim_values_.end())
    out = val->second;
  else {  // still trying to build it from (min/max) attributes
    fill<double>(key + "min", out.min());
    fill<double>(key + "max", out.max());
  }
  return out.validate();
}

template <>
const ParametersList& ParametersList::fill<Limits>(const std::string& key, Limits& value) const {
  if (has<Limits>(key)) {
    const auto& limits = get<Limits>(key);
    if (limits.hasMin())
      value.min() = limits.min();
    if (limits.hasMax())
      value.max() = limits.max();
    return *this;
  }
  fill<double>(key + "min", value.min());
  fill<double>(key + "max", value.max());
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
                                                           const ParticleProperties& default_value) const {
  if (has<ParametersList>(key)) {  // try to steer as a dictionary of particle properties
    const auto& plist = get<ParametersList>(key);
    if (plist.keys() == std::vector<std::string>{"pdgid"})
      return PDG::get()(plist.get<int>("pdgid"));
    return ParticleProperties(plist);
  }
  if (has<pdgid_t>(key)) {  // if not a dictionary of properties, retrieve from PDG runtime database
    CG_DEBUG("ParametersList") << "Retrieved physical properties for particle with PDG identifier '"
                               << get<pdgid_t>(key) << "' from PDG database.";
    return PDG::get()(get<pdgid_t>(key));
  }
  if (has<int>(key)) {  // if not a dictionary of properties, retrieve from PDG runtime database
    CG_DEBUG("ParametersList") << "Retrieved physical properties for particle with PDG identifier '" << get<int>(key)
                               << "' from PDG database.";
    return PDG::get()(get<int>(key));
  }
  CG_DEBUG("ParametersList") << "Failed to retrieve particle properties parameter with key=" << key << ".";
  return default_value;
}

/// Set a particle properties object value
template <>
ParametersList& ParametersList::set<ParticleProperties>(const std::string& key, const ParticleProperties& value) {
  PDG::get().define(value);
  return set(key, value.parameters());
}

template <>
std::vector<std::string> ParametersList::keysOf<ParticleProperties>() const {
  std::vector<std::string> steered_parameters_keys;
  for (const auto& key : keys())
    if (get<ParticleProperties>(key, ParticleProperties(-1)) != ParticleProperties(-1))
      steered_parameters_keys.emplace_back(key);
  return steered_parameters_keys;
}

namespace cepgen {
  std::ostream& operator<<(std::ostream& os, const ParametersList& parameters) {
    parameters.print(os);
    return os;
  }
}  // namespace cepgen

#undef IMPL_TYPE_SET
#undef IMPL_TYPE_GET
#undef IMPL_TYPE_ALL
