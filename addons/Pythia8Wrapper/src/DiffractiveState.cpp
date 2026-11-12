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
#include <Pythia8/Pythia.h>

#include "CepGenPythia8/DiffractiveState.h"
#include "CepGenPythia8/EventInterface.h"

using namespace cepgen::pythia8;

DiffractiveState::DiffractiveState(const Momentum& system,
                                   long beam_pdg_id,
                                   utils::RandomGenerator& rng,
                                   Pythia8::Pythia& pythia)
    : hadr::Hadroniser::DiffractiveState(system, beam_pdg_id), pythia_(&pythia) {
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
  // constituent masses, as used by the Lund string model
  valquark_mass_ = pythia_->particleData.m0(static_cast<int>(valquark_pdg_id_));
  diquark_mass_ = pythia_->particleData.m0(static_cast<int>(diquark_pdg_id_));

  // two-body decay momentum in the rest frame of the system
  const auto mass2 = system_mass_ * system_mass_;
  const auto sum = valquark_mass_ + diquark_mass_, diff = valquark_mass_ - diquark_mass_;
  const auto lambda = (mass2 - sum * sum) * (mass2 - diff * diff);
  if (lambda <= 0. || system_mass_ <= 0.) {
    CG_DEBUG("pythia8:DiffractiveState") << "Diffractive system mass unphysical. System invalid for fragmentation.";
    valid_ = false;
    return;
  }
  // isotropic string axis in dissociated system rest frame
  diquark_momentum_ =
      0.5 * std::sqrt(lambda) / system_mass_ * Momentum::randomThetaPhi(rng, {-M_PI, M_PI}, {0., 2. * M_PI});
}

bool DiffractiveState::hadronise(hadr::Hadroniser::Fragments& fragments, double tolerance) const {
  fragments.clear();

  // system fragmented at rest (q-qq string)
  constexpr int colour_tag = 101;
  const auto [quark_colour, quark_anticolour] =
      valquark_pdg_id_ > 0 ? std::make_pair(colour_tag, 0) : std::make_pair(0, colour_tag);
  const auto to_vec4 = [](const Momentum& p3, double mass) {
    return Pythia8::Vec4(p3.px(), p3.py(), p3.pz(), std::hypot(p3.p(), mass));
  };
  auto& event = pythia_->event;
  event.reset();
  event.append(static_cast<int>(valquark_pdg_id_),
               23,
               quark_colour,
               quark_anticolour,
               to_vec4(-diquark_momentum_, valquark_mass_),
               valquark_mass_);
  event.append(static_cast<int>(diquark_pdg_id_),
               23,
               quark_anticolour,
               quark_colour,
               to_vec4(diquark_momentum_, diquark_mass_),
               diquark_mass_);

  if (!pythia_->next())  // hadron-level standalone: string fragmentation and subsequent decays
    return false;

  for (int i = 3; i < event.size(); ++i) {  // skip the system entry and the two string ends
    const auto& particle = event[i];
    if (particle.statusAbs() < 81)  // only keep hadronisation (81-89) and decay (91-99) products
      continue;
    EventInterface::registerParticle(particle);
    fragments.emplace_back(hadr::Hadroniser::Fragment{
        i,
        particle.mother1(),  // string ends are not part of the fragments: attached to the system
        particle.id(),
        particle.isFinal(),
        Momentum::fromPxPyPzE(particle.px(), particle.py(), particle.pz(), particle.e())});
  }
  return momentumConserved(fragments, tolerance);
}
