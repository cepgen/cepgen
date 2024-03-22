/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2021-2024  Laurent Forthomme
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
#include "CepGen/Modules/ModuleFactoryImpl.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  template <typename T, typename I>
  ModuleFactory<T, I>::ModuleFactory(const std::string& descr) : description_(descr) {}

  template <typename T, typename I>
  std::unique_ptr<T> ModuleFactory<T, I>::build(const I& name, const ParametersList& params) const {
    if (name == I())
      throw CG_FATAL("ModuleFactory") << description_ << " cannot build a module with empty index/name!";
    auto plist = params;
    if (std::is_base_of<std::string, I>::value) {
      const auto extra_params = utils::split(utils::toString(name), '<');
      if (!extra_params.empty()) {
        plist.setName<std::string>(extra_params.at(0));
        if (extra_params.size() > 1)
          for (size_t i = 1; i < extra_params.size(); ++i)
            plist.feed(extra_params.at(i));
      }
    } else
      plist.setName<I>(name);
    return build(plist);
  }

  template <typename T, typename I>
  std::unique_ptr<T> ModuleFactory<T, I>::build(const ParametersList& params) const {
    if (!params.hasName<I>())
      throw CG_FATAL("ModuleFactory") << description_ << " failed to retrieve an indexing key "
                                      << "from parameters to build the module!\n"
                                      << "Parameters: " << params << ".\n"
                                      << "Registered modules: " << modules() << ".";
    const auto& idx = params.name<I>();
    if (map_.count(idx) == 0)
      throw CG_FATAL("ModuleFactory") << description_ << " failed to build a module with index/name \"" << idx
                                      << "\"!\nRegistered modules: " << modules() << ".";
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

  template <typename T, typename I>
  std::string ModuleFactory<T, I>::describe(const I& name) const {
    return describeParameters(name).description();
  }

  template <typename T, typename I>
  ParametersDescription ModuleFactory<T, I>::describeParameters(const ParametersList& params) const {
    if (!params.hasName<I>())
      throw CG_FATAL("ModuleFactory") << description_ << " failed to retrieve an indexing key "
                                      << "from parameters to describe the module!\n"
                                      << "Parameters: " << params << ".\n"
                                      << "Registered modules: " << modules() << ".";
    const auto& idx = params.name<I>();
    if (params_map_.count(idx) == 0)
      throw CG_FATAL("ModuleFactory") << "No parameters description were found for module index/name '" << idx << "'.\n"
                                      << "Registered modules: " << modules() << ".";
    return params_map_.at(idx).steer(params);
  }

  template <typename T, typename I>
  ParametersDescription ModuleFactory<T, I>::describeParameters(const I& name, const ParametersList& params) const {
    if (std::is_base_of<std::string, I>::value) {
      auto extra_params = utils::split(utils::toString(name), '<');
      auto* nm = reinterpret_cast<I*>(&extra_params[0]);
      if (params_map_.count(*nm) == 0)
        return ParametersDescription().setDescription("{module without description}").steer(params);
      auto descr = params_map_.at(*nm).steer(params);
      auto extra_params_obj = ParametersList();
      if (extra_params.size() > 1)
        for (size_t i = 1; i < extra_params.size(); ++i)
          extra_params_obj.feed(extra_params.at(i));
      return descr.steer(extra_params_obj);
    } else if (params_map_.count(name) == 0)
      return ParametersDescription().setDescription("{module without description}").steer(params);
    return params_map_.at(name).steer(params);
  }

  template <typename T, typename I>
  std::vector<I> ModuleFactory<T, I>::modules() const {
    std::vector<I> out;
    std::transform(map_.begin(), map_.end(), std::back_inserter(out), [](const auto& val) { return val.first; });
    std::sort(out.begin(), out.end());
    return out;
  }
}  // namespace cepgen
