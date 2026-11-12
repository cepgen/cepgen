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
#include <CepGen/Event/Event.h>
#include <CepGen/Physics/DiffractiveState.h>
#include <CepGen/Physics/PDG.h>
#include <Pythia8/Pythia.h>

#include "CepGenPythia8/EventInterface.h"

using namespace cepgen::pythia8;

namespace {
  /// Fraction of the decay channels left open for a given particle
  double openFraction(Pythia8::Pythia& pythia, const Pythia8::Particle& particle) {
    if (particle.isResonance())  // dynamic widths, computed by Pythia
      return pythia.particleData.resOpenFrac(particle.id());
    const auto& entry = particle.particleDataEntry();
    double fraction = 0.;
    for (int i = 0; i < entry.sizeChannels(); ++i)
      if (const auto& channel = entry.channel(i); channel.onMode() == 1 ||
                                                  (particle.id() > 0 && channel.onMode() == 2) ||
                                                  (particle.id() < 0 && channel.onMode() == 3))
        fraction += channel.bRatio();
    return fraction;
  }
}  // namespace

EventInterface::EventInterface(Event& event, Pythia8::Pythia& pythia) : cepgen_event_(event), pythia_(pythia) {}

void EventInterface::registerParticle(const Pythia8::Particle& particle) {
  if (PDG::get().has(particle.idAbs()))
    return;
  ParticleProperties properties;
  properties.pdgid = particle.idAbs();
  properties.name = properties.human_name = particle.name();
  properties.colours = particle.colType() == 0 ? 1 : std::abs(particle.colType()) == 1 ? 3 : 8;  // colour factor
  properties.mass = particle.m0();
  properties.width = particle.mWidth();
  if (const auto charge = particle.chargeType(); charge != 0)  // in units of e/3
    properties.charges = {charge, -charge};
  properties.fermion = particle.spinType() % 2 == 0;  // spinType = 2s+1
  PDG::get().define(properties);
}

bool EventInterface::fillEventBlock() {
  auto& event = pythia_.event;
  event.reset();
  cepgen_ids_.clear();

  // only the active particles are fed; Pythia 8 decays/hadronises anything left with a positive status.
  // coloured particles are sorted as (anti)triplets, octets, and antitriplets to build open strings
  std::vector<const Particle*> triplets, octets, antitriplets;
  for (const auto& particle : cepgen_event_(Particle::Role::CentralSystem)) {
    const auto status = particle.status();
    if (status != Particle::Status::FinalState && status != Particle::Status::Undecayed &&
        status != Particle::Status::Unfragmented)
      continue;  // already decayed/fragmented, or propagator
    const auto pdg_id = static_cast<int>(particle.integerPdgId());
    const auto colour_type = status == Particle::Status::Unfragmented ? pythia_.particleData.colType(pdg_id) : 0;
    if (colour_type == 1)
      triplets.emplace_back(&particle);
    else if (colour_type == -1)
      antitriplets.emplace_back(&particle);
    else if (colour_type == 2)
      octets.emplace_back(&particle);
    else {  // colour singlet
      const auto& mom = particle.momentum();
      cepgen_ids_[event.append(pdg_id, 23, 0, 0, mom.px(), mom.py(), mom.pz(), mom.energy(), mom.mass())] =
          particle.id();
    }
  }
  if (triplets.size() != antitriplets.size() || (triplets.empty() && octets.size() == 1)) {
    CG_WARNING("pythia8:EventInterface") << "Cannot build colour singlets out of " << triplets.size() << " triplet(s), "
                                         << antitriplets.size() << " antitriplet(s), and " << octets.size()
                                         << " octet(s) in the central system.";
    return false;
  }
  // gluons attached to first string, or form a closed loop if no triplets present
  int colour = 100;
  const auto append = [&](const Particle& particle, int col, int acol) {
    const auto& mom = particle.momentum();
    cepgen_ids_[event.append(
        static_cast<int>(particle.integerPdgId()), 23, col, acol, mom.px(), mom.py(), mom.pz(), mom.energy(), mom.mass())] =
        particle.id();
  };
  for (size_t i = 0; i < triplets.size(); ++i) {
    append(*triplets.at(i), ++colour, 0);
    if (i == 0)
      for (const auto* gluon : octets) {
        const auto anticolour = colour;
        append(*gluon, ++colour, anticolour);
      }
    append(*antitriplets.at(i), 0, colour);
  }
  if (triplets.empty() && !octets.empty()) {  // gluon loop
    const auto first_colour = colour + 1;
    for (size_t i = 0; i < octets.size(); ++i) {
      const auto anticolour = i == 0 ? first_colour + static_cast<int>(octets.size()) - 1 : colour;
      append(*octets.at(i), ++colour, anticolour);
    }
  }
  return true;
}

bool EventInterface::decay(double& weight) {
  if (!fillEventBlock())
    return false;
  const auto& event = pythia_.event;
  const auto old_particles_multiplicity = event.size();
  if (!pythia_.next()) {
    CG_WARNING("pythia8:EventInterface:decay") << "Pythia 8 decay/fragmentation of the central system failed.";
    return false;
  }

  for (const auto& [pythia_id, cepgen_id] : cepgen_ids_) {  // update fed particle status + branching fractions
    const auto& particle = event[pythia_id];
    if (particle.isFinal())
      continue;
    if (particle.statusAbs() >= 71 && particle.statusAbs() <= 79)  // hadronised in a string
      cepgen_event_[cepgen_id].setStatus(Particle::Status::Fragmented);
    else {  // decayed
      cepgen_event_[cepgen_id].setStatus(Particle::Status::Resonance);
      weight *= openFraction(pythia_, particle);
    }
  }

  for (int i = old_particles_multiplicity; i < event.size(); ++i) {  // append decay/hadronisation products to cent.sys.
    const auto& py_particle = event[i];
    registerParticle(py_particle);
    auto& particle = cepgen_event_.addParticle(Particle::Role::CentralSystem).get();
    particle.setIntegerPdgId(py_particle.id())
        .setStatus(py_particle.isFinal() ? Particle::Status::FinalState : Particle::Status::Resonance)
        .setMomentum(Momentum::fromPxPyPzE(py_particle.px(), py_particle.py(), py_particle.pz(), py_particle.e())
                         .setMass(py_particle.m()));
    const auto particle_id = particle.id();  // reference may be invalidated by further insertions
    for (const auto& mother_index : py_particle.motherList())
      if (const auto mother = cepgen_ids_.find(mother_index); mother != cepgen_ids_.end())
        cepgen_event_[particle_id].addMother(cepgen_event_[mother->second]);
    cepgen_ids_[i] = particle_id;
  }
  return true;
}

void EventInterface::update(Particle::Role role, const hadr::Fragments& fragments) {
  // particles appended to the vector the system lives in -> access through identifier
  const auto system_id = cepgen_event_.oneWithRole(role).id();
  const auto system_momentum = cepgen_event_[system_id].momentum();
  cepgen_event_[system_id].children().clear();  // events are recycled, with no reset of the system's filiation

  std::unordered_map<int, int> identifiers;  // Pythia -> CepGen
  for (const auto& fragment : fragments) {
    auto& particle = cepgen_event_.addParticle(role).get();
    particle.setIntegerPdgId(fragment.pdg_id)
        .setStatus(fragment.stable ? Particle::Status::FinalState : Particle::Status::Resonance)
        .setMomentum(Momentum(fragment.momentum).lorentzBoost(system_momentum), /*off_shell=*/true);
    const auto particle_id = particle.id();
    const auto mother = identifiers.find(fragment.mother);
    cepgen_event_[particle_id].addMother(cepgen_event_[mother != identifiers.end() ? mother->second : system_id]);
    identifiers[fragment.index] = particle_id;
  }
  cepgen_event_[system_id].setStatus(Particle::Status::Fragmented);  // after filiation (alters status)
}
