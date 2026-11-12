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
    /// Fill one end of the string in the Pythia 6 event record
    void fillParton(int index, long pdg_id, const Momentum& p3, double mass) {
      pyjets_.k[0][index] = 1;  // set to string flags by PYJOIN
      pyjets_.k[1][index] = static_cast<int>(pdg_id);
      pyjets_.k[2][index] = pyjets_.k[3][index] = pyjets_.k[4][index] = 0;
      pyjets_.p[0][index] = p3.px();
      pyjets_.p[1][index] = p3.py();
      pyjets_.p[2][index] = p3.pz();
      pyjets_.p[3][index] = std::hypot(p3.p(), mass);
      pyjets_.p[4][index] = mass;
      for (int i = 0; i < 5; ++i)
        pyjets_.v[i][index] = 0.;
    }
  }  // namespace
}  // namespace cepgen::pythia6

using namespace cepgen::pythia6;

DiffractiveState::DiffractiveState(const Momentum& system, long beam_pdg_id, utils::RandomGenerator& rng)
    : hadr::Hadroniser::DiffractiveState(system) {
  // valence flavour content, weighted by the quark electric charge squared and the diquark spin multiplicity
  const long sign = beam_pdg_id < 0 ? -1 : +1;
  if (const auto random_quarks = rng.uniformInt(0, 8); random_quarks < 1) {  // (d,uu1)
    valquark_pdg_id_ = sign * PDG::down;
    diquark_pdg_id_ = sign * PDG::uu1;
  } else if (random_quarks < 5) {  // (u,ud0)
    valquark_pdg_id_ = sign * PDG::up;
    diquark_pdg_id_ = sign * PDG::ud0;
  } else {  // (u,ud1)
    valquark_pdg_id_ = sign * PDG::up;
    diquark_pdg_id_ = sign * PDG::ud1;
  }
  checkPDGid(std::labs(valquark_pdg_id_));
  checkPDGid(std::labs(diquark_pdg_id_));
  valquark_mass_ = pymass(static_cast<int>(valquark_pdg_id_));
  diquark_mass_ = pymass(static_cast<int>(diquark_pdg_id_));

  // two-body decay momentum in the rest frame of the system
  const auto mass2 = system_mass_ * system_mass_;
  const auto sum = valquark_mass_ + diquark_mass_, diff = valquark_mass_ - diquark_mass_;
  const auto lambda = (mass2 - sum * sum) * (mass2 - diff * diff);
  if (lambda <= 0. || system_mass_ <= 0.) {
    CG_DEBUG("pythia6:DiffractiveState") << "Diffractive system mass unphysical. System invalid for fragmentation.";
    valid_ = false;
    return;
  }
  const auto p_star = 0.5 * std::sqrt(lambda) / system_mass_;

  // isotropic string axis in the rest frame of the dissociated system
  diquark_momentum_ = p_star * Momentum::randomThetaPhi(rng, {-M_PI, M_PI}, {0., 2. * M_PI});
}

bool DiffractiveState::hadronise(hadr::Hadroniser::Fragments& fragments, double tolerance) const {
  fragments.clear();

  // the system is fragmented at rest, as a single (anti)quark-(anti)diquark string
  pyjets_.n = 2;
  fillParton(0, valquark_pdg_id_, -diquark_momentum_, valquark_mass_);
  fillParton(1, diquark_pdg_id_, diquark_momentum_, diquark_mass_);
  pyjoin({1, 2});

  mstu(24) = 0;  // reset the latest error type
  pyexec();
  if (mstu(24) != 0)
    return false;

  for (int i = 2; i < pyjets_.n; ++i) {  // skip the two string ends
    const auto status = pyjets_.k[0][i];
    const bool stable = status == 1;
    if (!stable && status != 11)  // only keep stable and decayed particles
      continue;
    const auto pdg_id = pyjets_.k[1][i];
    if (const auto abs_pdg_id = std::abs(pdg_id); abs_pdg_id == 91 || abs_pdg_id == 92 || abs_pdg_id == 93)
      continue;  // skip the cluster/string/independent fragmentation system pseudo-particles
    checkPDGid(std::abs(pdg_id));
    fragments.emplace_back(hadr::Hadroniser::Fragment{
        i,
        pyjets_.k[2][i] - 1,
        pdg_id,
        stable,
        Momentum::fromPxPyPzE(pyjets_.p[0][i], pyjets_.p[1][i], pyjets_.p[2][i], pyjets_.p[3][i])});
  }
  return momentumConserved(fragments, tolerance);
}
