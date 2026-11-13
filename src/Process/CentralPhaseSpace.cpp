/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2026  Laurent Forthomme
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

#include <cmath>

#include "CepGen/Core/Exception.h"
#include "CepGen/Physics/Kinematics.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/CentralPhaseSpace.h"

using namespace cepgen;

CentralPhaseSpace::CentralPhaseSpace(const ParametersList& params) : NamedModule(params) {}

void CentralPhaseSpace::initialise(const Kinematics& kin, const spdgids_t& central_particles) {
  if (central_particles.size() != multiplicity())
    throw CG_FATAL("CentralPhaseSpace:initialise") << "Invalid central particles multiplicity: expecting "
                                                   << multiplicity() << ", got " << central_particles << ".";
  kin_ = &kin;
  particles_ = central_particles;
  masses_.clear();
  particle_cuts_.clear();
  const auto& per_id_cuts = kin.cuts().central_particles;
  for (const auto& pdg_id : particles_) {
    masses_.emplace_back(PDG::get().mass(pdg_id));
    // per-species cuts are steered for the unsigned PDG id (particle and antiparticle alike)
    if (const auto it = per_id_cuts.find(static_cast<pdgid_t>(std::abs(pdg_id))); it != per_id_cuts.end())
      particle_cuts_.emplace_back(it->second);
    else
      particle_cuts_.emplace_back();  // unconstrained
  }
}

const Kinematics& CentralPhaseSpace::kinematics() const {
  if (!kin_)
    throw CG_FATAL("CentralPhaseSpace:kinematics") << "Central phase space mapper was not yet initialised.";
  return *kin_;
}

ParametersDescription CentralPhaseSpace::description() {
  auto desc = NamedModule::description();
  desc.setDescription("Central system phase space mapper");
  return desc;
}
