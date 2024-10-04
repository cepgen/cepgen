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

#include "CepGen/Modules/PhaseSpaceGeneratorFactory.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Process/PhaseSpaceGenerator.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  PhaseSpaceGeneratorFactory& PhaseSpaceGeneratorFactory::get() {
    static PhaseSpaceGeneratorFactory instance;
    return instance;
  }

  std::unique_ptr<PhaseSpaceGenerator> PhaseSpaceGeneratorFactory::build(const ParametersList& params) const {
    if (const auto tokens = utils::split(params.name(), ':'); tokens.size() >= 2)
      return BasePhaseSpaceGeneratorFactory::build(
          ParametersList(params).setName(tokens.at(1)).set("partonsGenerator", tokens.at(0)));
    return BasePhaseSpaceGeneratorFactory::build(params);
  }
}  // namespace cepgen
