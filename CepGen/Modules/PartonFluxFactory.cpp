/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023  Laurent Forthomme
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
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  ParametersDescription PartonFluxFactory::describeParameters(const std::string& name,
                                                              const ParametersList& params) const {
    if (utils::contains(CollinearFluxFactory::get().modules(), name))
      return CollinearFluxFactory::get().describeParameters(name, params);
    if (utils::contains(KTFluxFactory::get().modules(), name))
      return KTFluxFactory::get().describeParameters(name, params);
    throw CG_FATAL("PartonFluxFactory:describeParameters")
        << "Failed to find a parton flux with name '" << name << "'.";
  }
}  // namespace cepgen
