/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2021-2025  Laurent Forthomme
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

#include <sstream>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersDescription.h"
#include "CepGen/Utils/Collections.h"
#include "CepGen/Utils/String.h"

using namespace cepgen;
using namespace std::string_literals;

ParametersDescription::ParametersDescription(const std::string& mod_key) {
  if (!mod_key.empty())
    setKey(mod_key);
}

ParametersDescription::ParametersDescription(const ParametersList& params) : ParametersList(params) {
  for (const auto& key : keys()) {
    if (obj_descr_.count(key) == 0)
      obj_descr_[key] = ParametersList::has<ParametersList>(key)
                            ? ParametersDescription(ParametersList::get<ParametersList>(key))
                            : ParametersDescription();
  }
  // avoid doubly-defined Limits/vector<double> situations
  for (const auto& limits_keys : keysOf<Limits>())
    if (utils::contains(keysOf<std::vector<double> >(),
                        limits_keys))  // hack to ensure vector<double> is dropped by set<Limits>(...)
      set(limits_keys, ParametersList::get<Limits>(limits_keys));
}

bool ParametersDescription::empty() const {
  // description is declared empty if
  return obj_descr_.empty() &&  // 1) it has no sub-description object, and
         mod_descr_.empty();    // 2) it is not itself describing anything
}

ParametersDescription& ParametersDescription::operator+=(const ParametersDescription& oth) {
  for (const auto& key : keysOf<std::vector<ParametersList> >())
    // particular case if one describes a set of key-indexed parameters list as a vector of parameters lists
    if (oth.parameters().has<ParametersList>(key)) {
      const auto& desc = get(key);
      erase(key);
      add(key, oth.get(key)).setDescription(desc.description()).allowedValues().append(desc.allowedValues());
    }
  obj_descr_.insert(oth.obj_descr_.begin(), oth.obj_descr_.end());
  obj_values_.append(oth.obj_values_);
  ParametersList::operator+=(oth);
  return *this;
}

bool ParametersDescription::has(const std::string& key) const { return obj_descr_.count(key) != 0; }

const ParametersDescription& ParametersDescription::get(const std::string& key) const {
  if (!has(key))
    throw CG_FATAL("ParametersDescription:get").log([&](auto& log) {
      log << "Failed to retrieve a parameters description member named \'" << key << "\'.\n"
          << "List of keys registered: ";
      std::string sep;
      for (const auto& [key, description] : obj_descr_)
        log << sep << "'" << key << "'", sep = ", ";
    });
  return obj_descr_.at(key);
}

ParametersDescription::Type ParametersDescription::type() const {
  if (is_vec_params_)
    return Type::ParametersVector;
  if (obj_descr_.empty())
    return Type::Value;
  if (const auto& mod_name = getNameString(); mod_name.empty())
    return Type::Parameters;
  return Type::Module;
}

std::string ParametersDescription::describe(size_t offset) const {
  static auto sep = [](size_t offset) -> std::string { return std::string(2 * offset, ' '); };
  const auto& mod_name = getNameString();
  const auto& parameters_description_type = type();
  const auto& keys = ParametersList::keys(false);
  std::ostringstream os;
  // write collection type (if collection)
  if (parameters_description_type == Type::Parameters)
    os << colourise("Parameters", utils::Colour::none, utils::Modifier::italic | utils::Modifier::underline)
       << " collection";
  else if (parameters_description_type == Type::Module)
    os << utils::boldify(mod_name) << " module";
  // write human-readable description (if exists)
  if (parameters_description_type != Type::ParametersVector && !mod_descr_.empty())
    os << " <- " << utils::colourise(mod_descr_, utils::Colour::blue, utils::Modifier::italic);
  if (keys.empty())  // no keys to this module ; can return
    return os.str();
  if (parameters_description_type == Type::ParametersVector) {
    os << parameters();
    return os.str();
  }
  if (parameters_description_type == Type::Module)
    os << " with parameters";
  os << ":";
  // write list of parameters (if it has some)
  for (const auto& key : keys) {
    os << "\n" << sep(offset + 1) << colourise(key, utils::Colour::none, utils::Modifier::underline) << " ";
    if (obj_descr_.count(key) == 0)
      continue;
    os << "= ";
    switch (const auto& parameter = obj_descr_.at(key); parameter.type()) {
      case Type::Value: {
        if (ParametersList::has<std::string>(key))
          os << "\"" << getString(key) << "\"";
        else
          os << getString(key);
        if (const auto& parameter_description = parameter.description(); !parameter_description.empty())
          os << " <- " << utils::colourise(parameter_description, utils::Colour::blue, utils::Modifier::italic);
        if (const auto& allowed_values = parameter.allowedValues(); !allowed_values.empty()) {
          os << " with allowed values:\n" << sep(offset + 2);
          std::string list_sep;
          for (const auto& [value, description] : allowed_values.allowed())
            os << list_sep << value << (!description.empty() ? " (" + description + ")" : ""), list_sep = ", ";
        }
      } break;
      case Type::ParametersVector: {
        os << colourise("Vector of parameters collections",
                        utils::Colour::none,
                        utils::Modifier::italic | utils::Modifier::underline);
        if (const auto& parameters_description = parameter.description(); !parameters_description.empty())
          os << " (" << colourise(parameters_description, utils::Colour::none, utils::Modifier::italic) << ")";
        if (const auto& parameters = ParametersList::get<ParametersList>(key); !parameters.empty())
          os << " with user-steered content: " << parameter.steer(parameters).describe(offset + 1);
        else
          os << " with expected content: " << parameter.describe(offset + 1);
      } break;
      default: {
        if (const auto& description = parameter.describe(offset + 1); !utils::trim(description).empty())
          os << description;
      } break;
    }
  }
  return os.str();
}

ParametersDescription& ParametersDescription::setDescription(const std::string& description) {
  mod_descr_ = description;
  return *this;
}

template <>
ParametersDescription& ParametersDescription::add<ParametersDescription>(const std::string& name,
                                                                         const ParametersDescription& desc) {
  obj_descr_[name] += desc;
  set<ParametersList>(name, desc.parameters());
  CG_DEBUG_LOOP("ParametersDescription:add").log([this, &name, &desc](auto& log) {
    log << "Added a new parameters collection \"" << name << "\" as: " << desc;
    const auto& mod_name = this->getNameString();
    if (!mod_name.empty())
      log << "\n"
          << "to the object with name: " << mod_name;
    log << ".";
  });
  return obj_descr_[name];
}

template <>
ParametersDescription& ParametersDescription::add<ParametersList>(const std::string& name, const ParametersList&) {
  throw CG_FATAL("ParametersDescription:add")
      << "Using a ParametersList object for the description of a collection of parameters is not allowed.\n"
      << "Please use a ParametersDescription object instead for the description of the '" << name << "' collection.";
}

ParametersDescription& ParametersDescription::addParametersDescriptionVector(
    const std::string& name,
    const ParametersDescription& description,
    const std::vector<ParametersList>& default_values) {
  obj_descr_[name] += description;
  obj_descr_[name].setParametersVector(true);
  auto& values = operator[]<std::vector<ParametersList> >(name);
  for (const auto& value : default_values)
    values.emplace_back(description.validate(value));
  CG_DEBUG_LOOP("ParametersDescription:addParametersDescriptionVector")
      .log([this, &name, &description, &default_values](auto& log) {
        log << "Added a new vector of parameters descriptions \"" << name << "\" as: " << description;
        if (const auto& mod_name = this->getNameString(); !mod_name.empty())
          log << "\nto the object with name: " << mod_name;
        log << ".\n";
        if (!default_values.empty())
          log << "It is now composed of " << default_values << ".";
      });
  return obj_descr_[name];
}

ParametersList& ParametersDescription::parameters() { return *this; }

const ParametersList& ParametersDescription::parameters() const { return *this; }

ParametersList ParametersDescription::validate(const ParametersList& user_params) const {
  auto plist = parameters();  // first copy the current parameters handled
  plist += user_params;
  for (const auto& key : keysOf<std::vector<ParametersList> >()) {
    if (user_params.has<std::vector<ParametersList> >(key)) {  // vector{ParametersList}
      plist.erase(key);
      for (const auto& pit : user_params.get<std::vector<ParametersList> >(key))
        plist.operator[]<std::vector<ParametersList> >(key).emplace_back(obj_descr_.at(key).parameters() + pit);
    } else if (user_params.has<ParametersList>(key)) {  // map{key -> ParametersList}
      plist.erase(key);
      const auto& pit = user_params.get<ParametersList>(key);
      for (const auto& kit : pit.keys())
        plist.operator[]<ParametersList>(key).operator[]<ParametersList>(kit) =
            obj_descr_.at(key).validate(pit.get<ParametersList>(kit));
    }
  }
  for (const auto& [key, parameters_description] : obj_descr_) {  // simple value
    if (!parameters_description.allowedValues().allAllowed()) {
      const auto validation_error =
          [](const auto& key, const auto& val, const ParametersDescription& desc) -> std::string {
        std::ostringstream os;
        os << "Invalid value for key '" << key << "'"
           << (!desc.description().empty() ? " ("s + desc.description() + ")" : "") << ".\n"
           << "Allowed values:\n";
        for (const auto& [value, description] : desc.allowedValues().allowed())
          os << " * " << value << (!description.empty() ? " (" + description + ")" : "") << "\n";
        os << "Steered value: '" << utils::toString(val) + "'.";
        return os.str();
      };
      if (plist.has<int>(key) && !parameters_description.allowedValues().validate(plist.get<int>(key)))
        throw CG_FATAL("ParametersDescription:validate")
            << validation_error(key, plist.get<int>(key), parameters_description);
      if (plist.has<std::string>(key) && !parameters_description.allowedValues().validate(plist.get<std::string>(key)))
        throw CG_FATAL("ParametersDescription:validate")
            << validation_error(key, plist.get<std::string>(key), parameters_description);
    }
  }
  CG_DEBUG_LOOP("ParametersDescription:validate") << "Validating user parameters:\n"
                                                  << "User-steered: " << user_params << ".\n"
                                                  << "Base/default: " << parameters() << ".\n"
                                                  << "-> Resulting: " << plist << ".\n";
  return plist;
}

ParametersDescription ParametersDescription::steer(const ParametersList& params) const {
  ParametersDescription parameters_description(*this);
  parameters_description += ParametersDescription(params);
  return parameters_description;
}

template <>
ParametersDescription& ParametersDescription::setKey<std::string>(const std::string& key) {
  mod_key_ = key;
  return *this;
}

ParametersDescription& ParametersDescription::allow(int val, const std::string& descr) {
  obj_values_.int_vals_[val] = descr;
  obj_values_.all_allowed_ = false;
  return *this;
}

ParametersDescription& ParametersDescription::allow(const std::string& val, const std::string& descr) {
  obj_values_.str_vals_[val] = descr;
  obj_values_.all_allowed_ = false;
  return *this;
}

bool ParametersDescription::ParameterValues::empty() const { return int_vals_.empty() && str_vals_.empty(); }

ParametersDescription::ParameterValues& ParametersDescription::ParameterValues::append(const ParameterValues& oth) {
  int_vals_.insert(oth.int_vals_.begin(), oth.int_vals_.end());
  str_vals_.insert(oth.str_vals_.begin(), oth.str_vals_.end());
  return *this;
}

std::map<std::string, std::string> ParametersDescription::ParameterValues::allowed() const {
  auto out = str_vals_;
  for (const auto& [value, description] : int_vals_)
    out[std::to_string(value)] = description;
  return out;
}

bool ParametersDescription::ParameterValues::validate(int val) const {
  if (allAllowed())
    return true;
  return int_vals_.count(val) > 0;
}

bool ParametersDescription::ParameterValues::validate(const std::string& val) const {
  if (allAllowed())
    return true;
  return str_vals_.count(val) > 0;
}

namespace cepgen {
  std::ostream& operator<<(std::ostream& os, const ParametersDescription& desc) { return os << desc.describe(); }

  std::ostream& operator<<(std::ostream& os, const ParametersDescription::Type& type) {
    switch (type) {
      case ParametersDescription::Type::Value:
        return os << "Value";
      case ParametersDescription::Type::Module:
        return os << "Module";
      case ParametersDescription::Type::Parameters:
        return os << "Parameters";
      case ParametersDescription::Type::ParametersVector:
        return os << "Parameters vector";
    }
    return os << "{invalid}";
  }

  std::ostream& operator<<(std::ostream& os, const ParametersDescription::ParameterValues& vals) {
    os << "Allowed values:";
    if (vals.allAllowed())
      return os << " any";
    if (vals.empty())
      return os << " none";
    if (!vals.int_vals_.empty())
      os << " integer{" << utils::repr(utils::keys(vals.int_vals_)) << "}";
    if (!vals.str_vals_.empty())
      os << " string{" << utils::repr(utils::keys(vals.str_vals_)) << "}";
    return os;
  }
}  // namespace cepgen
