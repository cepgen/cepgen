/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2025  Laurent Forthomme
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
#include "CepGen/Event/Event.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/Utils.h"
#include "CepGen/Utils/RandomGenerator.h"
#include "CepGen/Utils/String.h"
#include "CepGenPythia6/EventInterface.h"
#include "CepGenPythia6/Pythia6Interface.h"

using namespace cepgen::pythia6;

EventInterface::EventInterface(Event& event, const mode::Kinematics kinematics_mode, utils::RandomGenerator* rnd)
    : cepgen_event_(event), random_generator_(rnd) {
  if (kinematics_mode == mode::Kinematics::InelasticElastic || kinematics_mode == mode::Kinematics::InelasticInelastic)
    roles_.emplace_back(Particle::Role::OutgoingBeam1);
  if (kinematics_mode == mode::Kinematics::ElasticInelastic || kinematics_mode == mode::Kinematics::InelasticInelastic)
    roles_.emplace_back(Particle::Role::OutgoingBeam2);
}

void EventInterface::prepareHadronisation() {
  CG_DEBUG_LOOP("EventInterface:prepareHadronisation") << "Hadronisation preparation called.";

  for (auto role : roles_) {
    if (!cepgen_event_.hasRole(role))
      continue;
    const auto& part = cepgen_event_.oneWithRole(role);

    const auto [parton1, parton2] = pickPartonsContent();
    checkPDGid(parton1);
    checkPDGid(parton2);
    const double mq = pymass(parton1), mq2 = mq * mq;
    const double mdq = pymass(parton2), mdq2 = mdq * mdq;

    // choose random direction in MX frame
    const double phi = random_generator_->uniform(0., 2. * M_PI), theta = acos(random_generator_->uniform(-1., 1.));

    // compute momentum of decay particles from MX
    const double px2 = std::pow(utils::energyFromW(part.momentum().mass(), mdq2, mq2), 2) - mq2;
    if (px2 < 0.) {
      CG_WARNING("EventInterface:prepareHadronisation") << "Invalid remnants kinematics for " << role << ".";
      continue;
    }
    auto& beam_remnant = cepgen_event_[part.id()];

    // build 4-vectors and boost decay particles
    const double px = std::sqrt(px2);
    auto pdq = Momentum::fromPThetaPhiE(px, theta, phi, std::hypot(px, mdq)),
         pq = Momentum(-pdq).setEnergy(std::hypot(px, mq));

    // singlet
    cepgen_event_.addParticle(role)
        .get()
        .addMother(beam_remnant)
        .setPdgId(parton1, +1)
        .setStatus(Particle::Status::Unfragmented)
        .setMomentum(pq.lorentzBoost(beam_remnant.momentum()));

    // quark doublet
    cepgen_event_.addParticle(role)
        .get()
        .addMother(beam_remnant)
        .setPdgId(parton2, +1)
        .setStatus(Particle::Status::Unfragmented)
        .setMomentum(pdq.lorentzBoost(beam_remnant.momentum()));

    beam_remnant.setStatus(Particle::Status::Fragmented);
  }
}

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

      if (part.status() == Particle::Status::Unfragmented) {
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
      dbg << "Joining " << utils::s("particle", evt_string.size()) << " with " << cepgen_event_[evt_string[0]].role()
          << " role"
          << " in a same string";
      for (const auto& part_id : evt_string) {
        if (part_id != -1)
          dbg << utils::format("\n\t * %2d (pdgId=%4d)", part_id, pyjets_.k[1][part_id - 1]);
      }
    });
    pyjoin(evt_string);
  }
}

void EventInterface::run() {
  fillEventBlock();
  const auto old_particles_multiplicity = pyjets_.n;
  pyexec();
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
}

std::pair<short, short> EventInterface::pickPartonsContent() const {
  if (const auto random_quarks = random_generator_->uniformInt(0, 9); random_quarks < 1)
    return {PDG::down, 2203};  // (d,uu1)
  else if (random_quarks < 5)
    return {PDG::up, 2101};  // (u,ud0)
  return {PDG::up, 2103};    // (u,ud1)
}
