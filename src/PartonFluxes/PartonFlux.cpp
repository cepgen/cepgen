/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2025  Laurent Forthomme
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

#include "CepGen/PartonFluxes/PartonFlux.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/PDG.h"

using namespace cepgen;

PartonFlux::PartonFlux(const ParametersList& params)
    : NamedModule(params),
      alpha_over_pi_(constants::ALPHA_EM * M_1_PI),
      mp_(PDG::get().mass(PDG::proton)),
      mp2_(mp_ * mp_) {}

ParametersDescription PartonFlux::description() {
  auto desc = ParametersDescription();
  desc.setDescription("Unnamed parton flux evaluator");
  return desc;
}
