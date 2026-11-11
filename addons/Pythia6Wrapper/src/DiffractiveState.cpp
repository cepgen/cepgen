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

#include <CepGen/Core/Exception.h>
#include <CepGen/Physics/PDG.h>
#include <CepGen/Utils/RandomGenerator.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>

#include "CepGenPythia6/DiffractiveState.h"
#include "CepGenPythia6/Pythia6Interface.h"

namespace cepgen::pythia6 {
  namespace {
    constexpr long kDiquarkUU1 = 2203, kDiquarkUD0 = 2101, kDiquarkUD1 = 2103;

    /// Fill one end of the string in the Pythia 6 event record
    void fillParton(int index, long pdg_id, const std::array<double, 3>& p3, double mass) {
      pyjets_.k[0][index] = 1;  // set to string flags by PYJOIN
      pyjets_.k[1][index] = static_cast<int>(pdg_id);
      pyjets_.k[2][index] = pyjets_.k[3][index] = pyjets_.k[4][index] = 0;
      for (int i = 0; i < 3; ++i)
        pyjets_.p[i][index] = p3[i];
      pyjets_.p[3][index] = std::sqrt(p3[0] * p3[0] + p3[1] * p3[1] + p3[2] * p3[2] + mass * mass);
      pyjets_.p[4][index] = mass;
      for (int i = 0; i < 5; ++i)
        pyjets_.v[i][index] = 0.;
    }
  }  // namespace
}  // namespace cepgen::pythia6

using namespace cepgen::pythia6;

DiffractiveState DiffractiveState::build(const Momentum& system, long beam_pdg_id, utils::RandomGenerator& rng) {
  DiffractiveState diffractive_state;
  diffractive_state.system_mass_ = system.mass();

  // valence flavour content, weighted by the quark electric charge squared and the diquark spin multiplicity
  if (const auto random_quarks = rng.uniformInt(0, 8); random_quarks < 1)
    diffractive_state.struck_pdg_id_ = PDG::down, diffractive_state.spectator_pdg_id_ = kDiquarkUU1;  // (d,uu1)
  else if (random_quarks < 5)
    diffractive_state.struck_pdg_id_ = PDG::up, diffractive_state.spectator_pdg_id_ = kDiquarkUD0;  // (u,ud0)
  else
    diffractive_state.struck_pdg_id_ = PDG::up, diffractive_state.spectator_pdg_id_ = kDiquarkUD1;  // (u,ud1)
  if (beam_pdg_id < 0) {  // antibaryon content
    diffractive_state.struck_pdg_id_ *= -1;
    diffractive_state.spectator_pdg_id_ *= -1;
  }
  checkPDGid(std::labs(diffractive_state.struck_pdg_id_));
  checkPDGid(std::labs(diffractive_state.spectator_pdg_id_));
  diffractive_state.struck_mass_ = pymass(static_cast<int>(diffractive_state.struck_pdg_id_));
  diffractive_state.spectator_mass_ = pymass(static_cast<int>(diffractive_state.spectator_pdg_id_));

  // two-body decay momentum in the rest frame of the system
  const auto mass2 = diffractive_state.system_mass_ * diffractive_state.system_mass_;
  const auto sum = diffractive_state.struck_mass_ + diffractive_state.spectator_mass_,
             diff = diffractive_state.struck_mass_ - diffractive_state.spectator_mass_;
  const auto lambda = (mass2 - sum * sum) * (mass2 - diff * diff);
  if (lambda <= 0. || diffractive_state.system_mass_ <= 0.) {
    CG_DEBUG("pythia6:DiffractiveState") << "Diffractive system mass unphysical. System invalid for fragmentation.";
    diffractive_state.valid_ = false;
    return diffractive_state;
  }
  const auto p_star = 0.5 * std::sqrt(lambda) / diffractive_state.system_mass_;

  // isotropic string axis in the rest frame of the dissociated system
  const auto cos_theta = rng.uniform(-1., 1.), phi = rng.uniform(0., 2. * M_PI);
  const auto sin_theta = std::sqrt(1. - cos_theta * cos_theta);
  diffractive_state.spectator_momentum_ = {
      p_star * sin_theta * std::cos(phi), p_star * sin_theta * std::sin(phi), p_star * cos_theta};
  return diffractive_state;
}

bool DiffractiveState::hadronise(hadr::Hadroniser::Fragments& fragments, double tolerance) const {
  fragments.clear();

  // the system is fragmented at rest, as a single (anti)quark-(anti)diquark string
  const auto& spectator = spectator_momentum_;
  pyjets_.n = 2;
  fillParton(0, struck_pdg_id_, {-spectator[0], -spectator[1], -spectator[2]}, struck_mass_);
  fillParton(1, spectator_pdg_id_, spectator, spectator_mass_);
  pyjoin({1, 2});

  mstu(24) = 0;  // reset the latest error type
  pyexec();
  if (mstu(24) != 0)
    return false;

  double sum[4] = {0., 0., 0., 0.};
  for (int i = 2; i < pyjets_.n; ++i) {  // skip the two string ends
    const auto status = pyjets_.k[0][i];
    const bool stable = status == 1;
    if (!stable && status != 11)  // only keep stable and decayed particles
      continue;
    const auto pdg_id = pyjets_.k[1][i];
    if (const auto abs_pdg_id = std::abs(pdg_id); abs_pdg_id == 91 || abs_pdg_id == 92 || abs_pdg_id == 93)
      continue;  // skip the cluster/string/independent fragmentation system pseudo-particles
    checkPDGid(std::abs(pdg_id));
    hadr::Hadroniser::Fragment fragment;
    fragment.index = i;
    fragment.mother = pyjets_.k[2][i] - 1;
    fragment.pdg_id = pdg_id;
    fragment.stable = stable;
    fragment.momentum = Momentum::fromPxPyPzE(pyjets_.p[0][i], pyjets_.p[1][i], pyjets_.p[2][i], pyjets_.p[3][i]);
    fragments.emplace_back(fragment);
    if (stable)
      for (int j = 0; j < 4; ++j)
        sum[j] += pyjets_.p[j][i];
  }
  sum[3] -= system_mass_;
  const auto residual = std::max({std::fabs(sum[0]), std::fabs(sum[1]), std::fabs(sum[2]), std::fabs(sum[3])});
  if (residual > tolerance * system_mass_) {
    CG_WARNING("pythia6:hadronise") << "Momentum not conserved down to " << tolerance
                                    << " tolerance after hadronisation (residual=" << residual
                                    << " GeV for a system of mass " << system_mass_ << " GeV).";
    return false;
  }
  return true;
}
