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

#include <CepGen/Core/Exception.h>
#include <CepGen/Event/Event.h>
#include <CepGen/Utils/String.h>

#include <cstdlib>
#include <map>

#include "CepGenPythia6/EventInterface.h"
#include "CepGenPythia6/Pythia6Interface.h"

using namespace cepgen::pythia6;

namespace {
  bool isDissociatedBeam(const cepgen::Particle& particle) {
    return (particle.role() == cepgen::Particle::Role::OutgoingBeam1 ||
            particle.role() == cepgen::Particle::Role::OutgoingBeam2) &&
           particle.status() == cepgen::Particle::Status::Unfragmented;
  }
}  // namespace

EventInterface::EventInterface(Event& event) : cepgen_event_(event) {}

void EventInterface::fillEventBlock() {
  pyjets_.n = 0;         // reinitialise the event content
  evt_strings_.clear();  // reinitialise the string fragmentation variables

  for (const auto& role : cepgen_event_.roles()) {  // loop on roles
    string_t evt_string;
    for (const auto& part : cepgen_event_(role)) {
      const unsigned short i = part.id();
      pyjets_.p[0][i] = part.momentum().px();
      pyjets_.p[1][i] = part.momentum().py();
      pyjets_.p[2][i] = part.momentum().pz();
      pyjets_.p[3][i] = part.momentum().energy();
      pyjets_.p[4][i] = part.momentum().mass();
      try {
        pyjets_.k[0][i] = pythia6Status(static_cast<int>(part.status()));
      } catch (const std::out_of_range&) {
        cepgen_event_.dump();
        throw CG_FATAL("EventInterface") << "Failed to retrieve a Pythia 6 particle status translation for "
                                         << "CepGen status " << part.status() << "!";
      }
      pyjets_.k[1][i] = part.integerPdgId();
      const auto& moth = part.mothers();
      pyjets_.k[2][i] = moth.empty() ? 0                   // no mother
                                     : *moth.begin() + 1;  // mother
      const auto& daug = part.children();
      if (daug.empty())  // no children
        pyjets_.k[3][i] = pyjets_.k[4][i] = 0;
      else {
        pyjets_.k[3][i] = *daug.begin() + 1;   // first child
        pyjets_.k[4][i] = *daug.rbegin() + 1;  // last child
      }
      for (int j = 0; j < 5; ++j)
        pyjets_.v[j][i] = 0.;  // vertex position

      if (isDissociatedBeam(part))
        pyjets_.k[0][i] = 11;  // inactive entry, fragmented separately through a DiffractiveState
      else if (part.status() == Particle::Status::Unfragmented) {
        pyjets_.k[0][i] = 1;  // PYTHIA/JETSET workaround
        evt_string.emplace_back(part.id() + 1);
      } else if (part.status() == Particle::Status::Undecayed)
        pyjets_.k[0][i] = 2;  // intermediate resonance
      pyjets_.n++;
    }
    if (!evt_string.empty())  // at most one string per role
      evt_strings_.emplace_back(evt_string);
  }

  for (const auto& evt_string : evt_strings_) {  // loop over the strings to bind everything together
    if (evt_string.size() < 2)
      continue;

    CG_DEBUG_LOOP("EventInterface").log([&](auto& dbg) {
      dbg << "Joining " << utils::s("particle", evt_string.size()) << " with "
          << cepgen_event_[evt_string[0] - 1].role() << " role in a same string";
      for (const auto& particle_id : evt_string)
        dbg << utils::format("\n\t * %2d (pdgId=%4d)", particle_id, pyjets_.k[1][particle_id - 1]);
    });
    pyjoin(evt_string);
  }
}

bool EventInterface::decay() {
  fillEventBlock();
  const auto old_particles_multiplicity = pyjets_.n;
  mstu(24) = 0;  // reset the latest error type
  pyexec();
  if (mstu(24) != 0) {
    CG_WARNING("EventInterface:decay") << "Pythia 6 decay/fragmentation failed with error type " << mstu(24) << ".";
    return false;
  }
  //--- update the event
  for (int p = old_particles_multiplicity;  // filter the first particles already present in the event
       p < pyjets_.n;
       ++p) {
    checkPDGid(std::abs(pyjets_.k[1][p]));

    const auto moth_id = pyjets_.k[2][p] - 1;
    const auto role = pyjets_.k[2][p] != 0 ? cepgen_event_[moth_id].role()  // child particle inherits its mother's role
                                           : Particle::Role::UnknownRole;
    auto& particle =
        cepgen_event_.addParticle(role)
            .get()
            .setId(p)
            .setStatus(cepgenStatus(pyjets_.k[0][p]))
            .setIntegerPdgId(pyjets_.k[1][p])
            .setMomentum(
                Momentum(pyjets_.p[0][p], pyjets_.p[1][p], pyjets_.p[2][p], pyjets_.p[3][p]).setMass(pyjets_.p[4][p]));
    // define particle parentage
    auto& mother_particle = cepgen_event_[moth_id];
    if (role != Particle::Role::UnknownRole)
      mother_particle.setStatus(role == Particle::Role::CentralSystem ? Particle::Status::Resonance
                                                                      : Particle::Status::Fragmented);
    particle.addMother(mother_particle);
  }
  return true;
}

void EventInterface::update(Particle::Role role, const hadr::Fragments& fragments) {
  // particles are appended to the vector the system lives in, hence only access it through its identifier
  const auto system_id = cepgen_event_.oneWithRole(role).id();
  const auto system_momentum = cepgen_event_[system_id].momentum();
  cepgen_event_[system_id].children().clear();  // events are recycled, with no reset of the system's filiation

  std::map<int, int> identifiers;  // Pythia 6 event record position -> CepGen identifier
  for (const auto& fragment : fragments) {
    auto& particle = cepgen_event_.addParticle(role).get();
    particle.setIntegerPdgId(fragment.pdg_id)
        .setStatus(fragment.stable ? Particle::Status::FinalState : Particle::Status::Resonance)
        .setMomentum(Momentum(fragment.momentum).lorentzBoost(system_momentum), /*off_shell=*/true);
    const auto mother = identifiers.find(fragment.mother);
    particle.addMother(cepgen_event_[mother != identifiers.end() ? mother->second : system_id]);
    identifiers[fragment.index] = particle.id();
  }
  cepgen_event_[system_id].setStatus(Particle::Status::Fragmented);  // after the filiation, as it alters this status
}
