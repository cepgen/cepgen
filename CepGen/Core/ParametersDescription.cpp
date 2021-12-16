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

#include <sstream>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersDescription.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  ParametersDescription::ParametersDescription(const std::string& mod_name) {
    if (!mod_name.empty())
      setName(mod_name);
  }

  ParametersDescription::ParametersDescription(const ParametersList& params) : ParametersList(params) {
    for (const auto& key : ParametersList::keys()) {
      if (ParametersList::has<ParametersList>(key))
        obj_descr_[key] = ParametersDescription(ParametersList::get<ParametersList>(key));
      else
        obj_descr_[key] = ParametersDescription();
    }
    // avoid doubly-defined Limits/vector<double> situations
    for (const auto& klim : ParametersList::keysOf<Limits>())
      if (utils::contains(ParametersList::keysOf<std::vector<double> >(), klim))
        // hack to ensure vector<double> is dropped by set<Limits>(...)
        ParametersList::set<Limits>(klim, ParametersList::get<Limits>(klim));
  }

  bool ParametersDescription::empty() const { return obj_descr_.empty() && mod_descr_.empty(); }

  ParametersDescription& ParametersDescription::operator+=(const ParametersDescription& oth) {
    for (const auto& key : ParametersList::keysOf<std::vector<ParametersList> >())
      // particular case if one describes a set of key-indexed parameters list as a vector of parameters lists
      if (oth.has<ParametersList>(key)) {
        const auto desc = get(key).description();
        ParametersList::erase(key);
        add<ParametersDescription>(key, oth.get(key)).setDescription(desc);
      }
    obj_descr_.insert(oth.obj_descr_.begin(), oth.obj_descr_.end());
    ParametersList::operator+=(oth);
    return *this;
  }

  const ParametersDescription& ParametersDescription::get(const std::string& key) const {
    if (obj_descr_.count(key) == 0)
      throw CG_FATAL("ParametersDescription:get")
          << "Failed to retrieve a parameters description member named \'" << key << "\'!";
    return obj_descr_.at(key);
  }

  ParametersDescription::Type ParametersDescription::type() const {
    if (obj_descr_.empty())
      return Type::Value;
    const auto& mod_name = ParametersList::getString(ParametersList::MODULE_NAME);
    if (mod_name.empty())
      return Type::Parameters;
    return Type::Module;
  }

  std::string ParametersDescription::describe(size_t offset) const {
    static auto sep = [](size_t offset) -> std::string { return std::string(offset, '\t'); };
    const auto& mod_name = ParametersList::getString(ParametersList::MODULE_NAME);
    const auto& pdtype = type();
    const auto& keys = ParametersList::keys(false);
    std::ostringstream os;
    // write collection type (if collection)
    if (pdtype == Type::Parameters)
      os << utils::colourise("Parameters", utils::Colour::none, utils::Modifier::italic | utils::Modifier::underline)
         << " collection";
    else if (pdtype == Type::Module)
      os << utils::boldify(mod_name + " module");
    // write human-readable description (if exists)
    if (!mod_descr_.empty())
      os << " (" << utils::colourise(mod_descr_, utils::Colour::none, utils::Modifier::italic) << ")";
    if (!keys.empty())
      os << " with parameters:";
    // write list of parameters (if has some)
    for (const auto& key : keys) {
      os << "\n" << sep(offset + 1) << utils::colourise(key, utils::Colour::none, utils::Modifier::underline) << " ";
      if (obj_descr_.count(key) > 0) {
        os << "=";
        if (obj_descr_.at(key).type() == Type::Value)
          os << " " << ParametersList::getString(key);
        const auto& descr = obj_descr_.at(key).describe(offset + 1);
        if (!utils::trim(descr).empty())
          os << " " << descr;
      }
    }
    return os.str();
  }

  ParametersDescription& ParametersDescription::setDescription(const std::string& descr) {
    mod_descr_ = descr;
    return *this;
  }

  template <>
  ParametersDescription& ParametersDescription::add<ParametersDescription>(const std::string& name,
                                                                           const ParametersDescription& desc) {
    obj_descr_[name] = desc;
    ParametersList::set<ParametersList>(name, desc.parameters());
    CG_DEBUG("ParametersDescription:add").log([this, &name, &desc](auto& log) {
      log << "Added a new parameters collection \"" << name << "\" as: " << desc;
      const auto& mod_name = this->getString(ParametersList::MODULE_NAME);
      if (!mod_name.empty())
        log << "\nto the object with name: " << mod_name;
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

  ParametersDescription& ParametersDescription::addParametersDescriptionVector(const std::string& name,
                                                                               const ParametersDescription& desc) {
    obj_descr_[name] = desc;
    ParametersList::set<std::vector<ParametersList> >(name, {});
    CG_DEBUG("ParametersDescription:addParametersDescriptionVector").log([this, &name, &desc](auto& log) {
      log << "Added a new vector of parameters descriptions \"" << name << "\" as: " << desc;
      const auto& mod_name = this->getString(ParametersList::MODULE_NAME);
      if (!mod_name.empty())
        log << "\nto the object with name: " << mod_name;
      log << ".";
    });
    return obj_descr_[name];
  }

  ParametersList& ParametersDescription::parameters() { return *this; }

  const ParametersList& ParametersDescription::parameters() const { return *this; }

  ParametersList ParametersDescription::validate(const ParametersList& user_params) const {
    ParametersList plist = parameters();
    plist += user_params;
    for (const auto& key : keysOf<std::vector<ParametersList> >()) {
      if (user_params.has<std::vector<ParametersList> >(key)) {  // vector{ParametersList}
        plist.erase(key);
        for (const auto& pit : user_params.get<std::vector<ParametersList> >(key))
          plist.operator[]<std::vector<ParametersList> >(key).emplace_back(obj_descr_.at(key).parameters() + pit);
      } else if (user_params.has<ParametersList>(key)) {  // map{key -> ParametersList}
        plist.erase(key);
        const auto& pit = user_params.get<ParametersList>(key);
        for (const auto& kit : pit.keys()) {
          plist.operator[]<ParametersList>(key).set<ParametersList>(
              kit, obj_descr_.at(key).parameters() + pit.get<ParametersList>(kit));
        }
      }
    }
    CG_DEBUG("ParametersDescription:validate") << "Validating user-defined parameters:\n"
                                               << user_params << ".\nBase parameters:\n"
                                               << parameters() << ".\nResult:\n"
                                               << plist << ".";
    return plist;
  }

  ParametersDescription ParametersDescription::steer(const ParametersList& params) const {
    ParametersDescription pdesc(*this);
    pdesc += ParametersDescription(params);
    return pdesc;
  }

  template <>
  ParametersDescription& ParametersDescription::setName<std::string>(const std::string& name) {
    mod_name_ = name;
    add<std::string>(ParametersList::MODULE_NAME, name);
    return *this;
  }

  std::ostream& operator<<(std::ostream& os, const ParametersDescription& desc) { return os << desc.describe(); }

  std::ostream& operator<<(std::ostream& os, const ParametersDescription::Type& type) {
    switch (type) {
      case ParametersDescription::Type::Value:
        return os << "Value";
      case ParametersDescription::Type::Module:
        return os << "Module";
      case ParametersDescription::Type::Parameters:
        return os << "Parameters";
    }
    return os << "{invalid}";
  }
}  // namespace cepgen
