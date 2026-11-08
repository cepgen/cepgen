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

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/RunParameters.h"
#include "CepGen/Modules/ModuleFactory.h"
#include "CepGen/Utils/String.h"

using namespace cepgen;

template <typename T>
ModuleFactory<T>::ModuleFactory(const std::string& description) : description_(description) {}

template <typename T>
std::unique_ptr<T> ModuleFactory<T>::build(const ParametersList& params) const {
  if (!params.hasName())
    throw CG_FATAL("ModuleFactory") << description_ << " failed to retrieve an indexing key "
                                    << "from parameters to build the module!\n"
                                    << "Parameters: " << params << ".\n"
                                    << "Registered modules: " << modules() << ".";
  const auto& idx = params.name();
  if (map_.count(idx) == 0)
    throw CG_FATAL("ModuleFactory") << description_ << " failed to build a module with name '" << idx << "'.\n"
                                    << "Registered modules: " << modules() << ".";
  ParametersList plist(describeParameters(idx).validate(params));
  CG_DEBUG("ModuleFactory").log([&](auto& log) {
    log << description_ << " will build a module ";
    if (plist.empty())
      log << "without parameters.";
    else
      log << "with parameters:\n" << plist << ".";
  });
  return map_.at(idx)(plist);
}

template <typename T>
std::unique_ptr<T> ModuleFactory<T>::build(const std::string& name, const ParametersList& params) const {
  if (name.empty())
    throw CG_FATAL("ModuleFactory") << description_ << " cannot build a module with empty index/name!";
  auto plist = params;
  if (const auto extra_params = utils::split(name, '<'); !extra_params.empty()) {
    plist.setName(extra_params.at(0));
    if (extra_params.size() > 1)
      for (size_t i = 1; i < extra_params.size(); ++i)
        plist.feed(extra_params.at(i));
  }
  return build(plist);
}

template <typename T>
std::unique_ptr<T> ModuleFactory<T>::build(int index, const ParametersList& params) const {
  if (indices_.count(index) > 0)
    return build(indices_.at(index), params);
  const auto& mod_names = modules();
  if (const auto str_index = std::to_string(index);
      std::find(mod_names.begin(), mod_names.end(), str_index) != mod_names.end())
    return build(str_index, params);
  throw CG_FATAL("ModuleFactory") << description_ << " failed to build a module with index '" << index << "'. \n"
                                  << "Registered indices: " << indices_ << ".";
}

template <typename T>
std::string ModuleFactory<T>::describe(const std::string& name) const {
  return describeParameters(name).description();
}

template <typename T>
ParametersDescription ModuleFactory<T>::describeParameters(const ParametersList& parameters) const {
  if (!parameters.hasName())
    throw CG_FATAL("ModuleFactory") << description_ << " failed to retrieve an indexing key "
                                    << "from parameters to describe the module!\n"
                                    << "Parameters: " << parameters << ".\n"
                                    << "Registered modules: " << modules() << ".";
  const auto& idx = parameters.name();
  if (params_map_.count(idx) == 0)
    throw CG_FATAL("ModuleFactory") << "No parameters description were found for module name '" << idx << "'.\n"
                                    << "Registered modules: " << modules() << ".";
  return params_map_.at(idx).steer(parameters);
}

template <typename T>
ParametersDescription ModuleFactory<T>::describeParameters(const std::string& name,
                                                           const ParametersList& params) const {
  const auto extra_params = utils::split(name, '<');
  const auto mod_name = extra_params.at(0);
  if (params_map_.count(mod_name) == 0)
    return ParametersDescription().setName(mod_name).setDescription("{module without description}").steer(params);
  auto description = params_map_.at(mod_name).steer(params);
  auto extra_params_obj = ParametersList();
  if (extra_params.size() > 1)
    for (size_t i = 1; i < extra_params.size(); ++i)
      extra_params_obj.feed(extra_params.at(i));
  return description.steer(extra_params_obj);
}

template <typename T>
ParametersDescription ModuleFactory<T>::describeParameters(int index, const ParametersList& params) const {
  if (indices_.count(index) > 0)
    return describeParameters(indices_.at(index), params);
  const auto& mod_names = modules();
  if (const auto str_index = std::to_string(index);
      std::find(mod_names.begin(), mod_names.end(), str_index) != mod_names.end())
    return describeParameters(str_index, params);
  throw CG_FATAL("ModuleFactory") << "No parameters description were found for module index '" << index << "'.\n"
                                  << "Registered modules: " << indices_ << ".";
}

template <typename T>
std::vector<std::string> ModuleFactory<T>::modules() const {
  std::vector<std::string> out;
  std::transform(map_.begin(), map_.end(), std::back_inserter(out), [](const auto& val) { return val.first; });
  std::sort(out.begin(), out.end());
  return out;
}
#include "CepGen/Modules/ModuleFactoryImpl.h"
