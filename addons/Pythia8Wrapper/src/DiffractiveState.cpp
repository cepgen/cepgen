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
#include <CepGen/Event/Particle.h>
#include <CepGen/Physics/PDG.h>
#include <CepGen/Utils/RandomGenerator.h>
#include <Pythia8/Pythia.h>

#include "CepGenPythia8/DiffractiveState.h"
#include "CepGenPythia8/EventInterface.h"

using namespace cepgen::pythia8;

DiffractiveState::DiffractiveState(utils::RandomGenerator& rng,
                                   const Momentum& system,
                                   const Particle& beam,
                                   Pythia8::Pythia& pythia,
                                   double primordial_kt)
    : hadr::Hadroniser::DiffractiveState(
          rng, beam.integerPdgId(), {{PDG::down, PDG::uu1, 1.}, {PDG::up, PDG::ud0, 4.}, {PDG::up, PDG::ud1, 4.}}),
      pythia_(&pythia) {
  initialise(system, beam.momentum(), primordial_kt);
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

double DiffractiveState::massFromPdgId(long pdg_id) const {
  if (pythia_)
    return pythia_->particleData.m0(static_cast<int>(pdg_id));
  throw CG_FATAL("pythia6:DiffractiveState:massFromPdgId")
      << "Trying to retrieve particle mass while no Pythia 8 object was initialised.";
}
