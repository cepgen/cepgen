/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2026  Laurent Forthomme
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

#include "CepGen/Modules/PartonsPhaseSpaceGeneratorFactory.h"
#include "CepGen/Process/PartonSampler.h"
#include "CepGen/Process/PartonsPhaseSpaceGenerator.h"

using namespace cepgen;

/// Collinear two-parton phase space generator: one collinear sampler per incoming beam
class PartonsCollinearPhaseSpaceGenerator final : public PartonsPhaseSpaceGenerator {
public:
  explicit PartonsCollinearPhaseSpaceGenerator(const ParametersList& params)
      : PartonsPhaseSpaceGenerator(params,
                                   std::make_unique<CollinearPartonSampler>(params),
                                   std::make_unique<CollinearPartonSampler>(params)) {}

  static ParametersDescription description() {
    auto desc = PartonsPhaseSpaceGenerator::description();
    desc.setDescription("Collinear phase space mapper");
    return desc;
  }
};
REGISTER_PARTONS_PHASE_SPACE_GENERATOR("coll", PartonsCollinearPhaseSpaceGenerator);
