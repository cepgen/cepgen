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

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/ModuleFactory.h"

// collection of handled objects
#include "CepGen/Cards/Handler.h"
#include "CepGen/Core/EventModifier.h"
#include "CepGen/Core/ExportModule.h"
#include "CepGen/Event/Event.h"
#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Physics/Coupling.h"
#include "CepGen/Processes/Process.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/StructureFunctions/SigmaRatio.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Functional.h"

namespace cepgen {
  template <typename T, typename I>
  ModuleFactory<T, I>& ModuleFactory<T, I>::get() {
    static ModuleFactory<T, I> instance;
    return instance;
  }

  /*template <typename T, typename I>
  template <typename U>
  void ModuleFactory<T, I>::registerModule(const I& name, const ParametersList& def_params) {
    static_assert(std::is_base_of<T, U>::value,
                  "\n\n  *** Failed to register an object with improper inheritance into the factory. ***\n");
    if (has(name))
      throw CG_FATAL("ModuleFactory") << description_ << " detected a duplicate module registration for index/name \""
                                      << name << "\"!";
    map_[name] = &build<U>;
    descr_map_[name] = U::description();
    params_map_[name] = !def_params.empty() ? ParametersDescription(def_params) : U::description();
    params_map_[name].setName(name);
  }*/

  template <typename T, typename I>
  std::unique_ptr<T> ModuleFactory<T, I>::build(const I& name, const ParametersList& params) const {
    if (name == I())
      throw CG_FATAL("ModuleFactory") << description_ << " cannot build a module with empty index/name!";
    return build(ParametersList(params).setName<I>(name));
  }

  template <typename T, typename I>
  std::unique_ptr<T> ModuleFactory<T, I>::build(const ParametersList& params) const {
    if (!params.has<I>(ParametersList::MODULE_NAME))
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
  ParametersDescription ModuleFactory<T, I>::describeParameters(const I& name, const ParametersList& params) const {
    if (params_map_.count(name) == 0)
      return ParametersDescription().setDescription("{module without description}").steer(params);
    return params_map_.at(name).steer(params);
  }

  template <typename T, typename I>
  std::vector<I> ModuleFactory<T, I>::modules() const {
    std::vector<I> out;
    std::transform(map_.begin(), map_.end(), std::back_inserter(out), [](auto& val) { return val.first; });
    std::sort(out.begin(), out.end());
    return out;
  }

  template class ModuleFactory<card::Handler, std::string>;
  template class ModuleFactory<Coupling, std::string>;
  template class ModuleFactory<utils::Drawer, std::string>;
  template class ModuleFactory<EventModifier, std::string>;
  template class ModuleFactory<formfac::Parameterisation, std::string>;
  template class ModuleFactory<Integrator, std::string>;
  template class ModuleFactory<io::ExportModule, std::string>;
  template class ModuleFactory<proc::Process, std::string>;
  template class ModuleFactory<sigrat::Parameterisation, int>;
  template class ModuleFactory<strfun::Parameterisation, int>;
  template class ModuleFactory<utils::Functional, std::string>;
}  // namespace cepgen
