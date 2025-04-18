/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2017-2025  Laurent Forthomme
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

#include <algorithm>
#include <cmath>

#include "CepGen/Event/Event.h"
#include "CepGen/Physics/Cuts.h"
#include "CepGen/Utils/Message.h"
#include "CepGen/Utils/String.h"

using namespace cepgen;
using namespace cepgen::cuts;

Central::Central() : SteeredObject(ParametersList{}) {}

Central::Central(const ParametersList& params) : SteeredObject(params) {
  (*this)
      .add("pt", pt_single)
      .add("eta", eta_single)
      .add("phi", phi_single)
      .add("rapidity", rapidity_single)
      .add("energy", energy_single)
      .add("mass", mass_single)
      .add("ptsum", pt_sum)
      .add("etasum", eta_sum)
      .add("energysum", energy_sum)
      .add("invmass", mass_sum)
      .add("ptdiff", pt_diff)
      .add("dphi", phi_diff)
      .add("rapiditydiff", rapidity_diff);
  if (params.has<Limits>("phiptdiff")) {
    CG_WARNING("Central") << "\"phiptdiff\" parameter is deprecated! "
                          << "Please use \"phidiff\" instead.";
    params.fill("phiptdiff", phi_diff);  // legacy
  }
}

bool Central::contain(const Particles& parts, const Event*) const {
  Momentum mom_sum;
  for (const auto& part : parts) {
    const auto& mom = part.momentum();
    if (!pt_single.contains(mom.pt()) || !eta_single.contains(mom.eta()) || !rapidity_single.contains(mom.rapidity()) ||
        !energy_single.contains(mom.energy()) || !mass_single.contains(mom.mass()))
      return false;
    mom_sum += mom;
  }
  if (!pt_sum.contains(mom_sum.pt()) || !eta_sum.contains(mom_sum.eta()) || !energy_sum.contains(mom_sum.energy()) ||
      !mass_sum.contains(mom_sum.mass()))
    return false;
  if (parts.size() > 1) {  // look at correlations
    const auto &mom1 = parts.at(0).momentum(), &mom2 = parts.at(1).momentum();
    if (!pt_diff.contains(std::fabs(mom1.pt() - mom2.pt())) || !phi_diff.contains(mom1.deltaPhi(mom2)) ||
        !rapidity_diff.contains(std::fabs(mom1.rapidity() - mom2.rapidity())))
      return false;
  }
  return true;
}

ParametersDescription Central::description() {
  auto desc = ParametersDescription();
  desc.add("pt", Limits{0.}).setDescription("Single particle pt (GeV/c)");
  desc.add("eta", Limits{}).setDescription("Single particle eta");
  desc.add("phi", Limits{0., 2. * M_PI}).setDescription("Single particle azimuthal angle");
  desc.add("rapidity", Limits{}).setDescription("Single particle rapidity");
  desc.add("energy", Limits{}).setDescription("Single particle energy (GeV)");
  desc.add("mass", Limits{}).setDescription("Single particle mass (GeV/c^2)");
  desc.add("ptsum", Limits{}).setDescription("System pt (GeV/c)");
  desc.add("etasum", Limits{}).setDescription("System eta");
  desc.add("energysum", Limits{}).setDescription("System energy (GeV)");
  desc.add("invmass", Limits{}).setDescription("System mass (GeV/c^2)");
  desc.add("ptdiff", Limits{}).setDescription("System D(pt) (GeV/c)");
  desc.add("dphi", Limits{}).setDescription("System D(phi) (rad)");
  desc.add("rapiditydiff", Limits{}).setDescription("System D(Y)");
  return desc;
}

//------------------------------------------------------------------

Initial::Initial(const ParametersList& params) : SteeredObject(params) {
  (*this).add("q2", q2).add("qt", qt).add("phi", phi);
}

void Initial::setParameters(const ParametersList& params) {
  if (params.empty())
    return;
  SteeredObject::setParameters(params);
  if (const auto q2lim = params.get<Limits>("q2"); q2lim.valid())  // symmetric Q^2 cut specified
    q2 = {q2lim, q2lim};
  for (auto& q2lim : q2)
    if (q2lim.max() <= 0.) {
      CG_WARNING("Initial:setParameters") << "Maximum parton virtuality (" << q2lim << ") is invalid. "
                                          << "It is now set to " << 1.e4 << " GeV^2.";
      q2lim.max() = 1.e4;
    }
}

bool Initial::contain(const Particles& parts, const Event*) const {
  for (const auto& part : parts) {
    const auto& mom = part.momentum();
    if (!qt.contains(mom.pt()))
      return false;
  }
  if (parts.size() == 2) {
    for (size_t i = 0; i < 2; ++i)
      if (!q2.at(i).contains(parts.at(i).momentum().mass2()))
        return false;
    if (phi.valid() && !phi.contains(parts.at(0).momentum().deltaPhi(parts.at(1).momentum())))
      return false;
  }
  return true;
}

ParametersDescription Initial::description() {
  auto desc = ParametersDescription();
  desc.add("q2", std::vector<Limits>(2, {0., 1.e5})).setDescription("Parton virtuality(ies) (GeV^2)");
  desc.add("qt", Limits{}).setDescription("Transverse virtuality (GeV)");
  desc.add("phi", Limits{}).setDescription("Partons D(phi) (rad)");
  return desc;
}

//------------------------------------------------------------------

Remnants::Remnants(const ParametersList& params) : SteeredObject(params) {
  (*this).add("mx", mx).add("yj", yj).add("xi", xi);
}

void Remnants::setParameters(const ParametersList& params) {
  if (params.empty())
    return;
  SteeredObject::setParameters(params);
  if (mx.min() < MX_MIN) {
    CG_WARNING("CutsList:setParameters") << "Minimum diffractive mass range (" << mx << ") is invalid. "
                                         << "It is now set to " << MX_MIN << " GeV/c^2.";
    mx.min() = MX_MIN;
  }
}

bool Remnants::contain(const Particles& parts, const Event* evt) const {
  for (const auto& part : parts) {
    if (part.status() != Particle::Status::FinalState)
      continue;
    if (evt && xi.valid() && !xi.contains(1. - part.momentum().pz() / (*evt)(*part.mothers().begin()).momentum().pz()))
      return false;
    if (!yj.contains(std::fabs(part.momentum().rapidity())))
      return false;
  }
  return true;
}

ParametersDescription Remnants::description() {
  auto desc = ParametersDescription();
  desc.add("mx", Limits{MX_MIN, 1.e3}).setDescription("Diffractive mass (GeV/c^2)");
  desc.add("yj", Limits{}).setDescription("Diffractive jet rapidity");
  desc.add("xi", Limits{}).setDescription("Longit.fract.mom. loss (\"xi\")");
  return desc;
}
