/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2026  Laurent Forthomme
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

#include "CepGen/Core/ParametersList.h"
#include "CepGen/Physics/Hadroniser.h"
#include "CepGen/Utils/Message.h"

using namespace cepgen;
using namespace cepgen::hadr;

Hadroniser::Hadroniser(const ParametersList& plist)
    : EventModifier(plist), fragment_remnants_(steer<bool>("remnantsFragmentation")) {}

ParametersDescription Hadroniser::description() {
  auto desc = EventModifier::description();
  desc.add("remnantsFragmentation", true).setDescription("Apply the fragmentation algorithm to proton remnants");
  return desc;
}

Hadroniser::DiffractiveState::DiffractiveState(const Momentum& system_momentum)
    : system_mass_(system_momentum.mass()) {}

bool Hadroniser::DiffractiveState::momentumConserved(const Fragments& fragments, double tolerance) const {
  Momentum sum{0., 0., 0., system_mass_};
  for (const auto& fragment : fragments)
    if (fragment.stable)
      sum -= fragment.momentum;
  if (const auto residual = sum.p(); residual > tolerance * system_mass_) {
    CG_WARNING("DiffractiveState:momentumConserved")
        << "Momentum not conserved down to " << tolerance << " tolerance after hadronisation (residual=" << residual
        << " GeV for a system of mass " << system_mass_ << " GeV).";
    return false;
  }
  return true;
}
