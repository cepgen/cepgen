/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2024  Laurent Forthomme
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

#include "CepGen/CollinearFluxes/CollinearFlux.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/KTFluxes/KTFlux.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Utils/Collections.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  ParametersDescription PartonFluxFactory::describeParameters(const std::string& name,
                                                              const ParametersList& params) const {
    if (name.empty())
      throw CG_FATAL("PartonFluxFactory:describeParameters") << "No name given to describe parton flux modelling.";
    if (utils::contains(CollinearFluxFactory::get().modules(), name))
      return CollinearFluxFactory::get().describeParameters(name, params);
    if (utils::contains(KTFluxFactory::get().modules(), name))
      return KTFluxFactory::get().describeParameters(name, params);
    return ParametersDescription().setName(name);
  }

  bool PartonFluxFactory::elastic(const ParametersList& params) const {
    const auto& name = params.name();
    if (name.empty())
      throw CG_FATAL("PartonFluxFactory:elastic") << "No name given to get parton flux modelling elasticity.";
    if (utils::contains(CollinearFluxFactory::get().modules(), name))
      return !CollinearFluxFactory::get().build(name, params)->fragmenting();
    if (utils::contains(KTFluxFactory::get().modules(), name))
      return !KTFluxFactory::get().build(name, params)->fragmenting();
    throw CG_FATAL("PartonFluxFactory:elastic") << "Failed to find a parton flux with name '" << name << "'.";
  }

  int PartonFluxFactory::partonPdgId(const ParametersList& params) const {
    const auto& name = params.name();
    if (name.empty())
      throw CG_FATAL("PartonFluxFactory:partonPdgId") << "No name given to get parton flux modelling PDG id.";
    if (utils::contains(CollinearFluxFactory::get().modules(), name))
      return CollinearFluxFactory::get().build(name, params)->partonPdgId();
    if (utils::contains(KTFluxFactory::get().modules(), name))
      return KTFluxFactory::get().build(name, params)->partonPdgId();
    throw CG_FATAL("PartonFluxFactory:partonPdgId") << "Failed to find a parton flux with name '" << name << "'.";
  }
}  // namespace cepgen
