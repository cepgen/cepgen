/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2024  Laurent Forthomme
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

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Event/Particle.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/Utils.h"
#include "CepGen/Utils/RandomGenerator.h"
#include "CepGen/Utils/String.h"
#include "CepGenPythia6/EventInterface.h"
#include "CepGenPythia6/Pythia6Interface.h"

namespace cepgen::pythia6 {
  EventInterface::EventInterface(Event& event, const mode::Kinematics kin_mode, utils::RandomGenerator* rnd)
      : evt_(event), rnd_(rnd) {
    if (kin_mode == mode::Kinematics::InelasticElastic || kin_mode == mode::Kinematics::InelasticInelastic)
      roles_.emplace_back(Particle::Role::OutgoingBeam1);
    if (kin_mode == mode::Kinematics::ElasticInelastic || kin_mode == mode::Kinematics::InelasticInelastic)
      roles_.emplace_back(Particle::Role::OutgoingBeam2);
  }

  void EventInterface::prepareHadronisation() {
    CG_DEBUG_LOOP("EventInterface:prepareHadronisation") << "Hadronisation preparation called.";

    for (auto role : roles_) {
      if (!evt_.hasRole(role))
        continue;
      const auto& part = evt_.oneWithRole(role);

      const auto partons = pickPartonsContent();
      checkPDGid(partons.first);
      checkPDGid(partons.second);
      const double mq = pymass(partons.first), mq2 = mq * mq;
      const double mdq = pymass(partons.second), mdq2 = mdq * mdq;

      // choose random direction in MX frame
      const double phi = rnd_->uniform(0., 2. * M_PI), theta = acos(rnd_->uniform(-1., 1.));

      // compute momentum of decay particles from MX
      const double px2 = std::pow(utils::energyFromW(part.momentum().mass(), mdq2, mq2), 2) - mq2;
      if (px2 < 0.) {
        CG_WARNING("EventInterface:prepareHadronisation") << "Invalid remnants kinematics for " << role << ".";
        continue;
      }
      const double px = std::sqrt(px2);

      //--- build 4-vectors and boost decay particles
      auto pdq = Momentum::fromPThetaPhiE(px, theta, phi, std::hypot(px, mdq));
      auto pq = -pdq;
      pq.setEnergy(std::hypot(px, mq));

      const auto part_id = part.id();

      // singlet
      auto& quark = evt_.addParticle(role).get();
      quark.addMother(evt_[part_id]);
      quark.setPdgId(partons.first, +1);
      quark.setStatus(Particle::Status::Unfragmented);
      quark.setMomentum(pq.lorentzBoost(evt_[part_id].momentum()));

      // quark doublet
      auto& diquark = evt_.addParticle(role).get();
      diquark.addMother(evt_[part_id]);
      diquark.setPdgId(partons.second, +1);
      diquark.setStatus(Particle::Status::Unfragmented);
      diquark.setMomentum(pdq.lorentzBoost(evt_[part_id].momentum()));

      evt_[part_id].setStatus(Particle::Status::Fragmented);
    }
  }

  void EventInterface::fillEventBlock() {
    pyjets_.n = 0;         // reinitialise the event content
    evt_strings_.clear();  // reinitialise the string fragmentation variables

    for (const auto& role : evt_.roles()) {  // loop on roles
      string_t evt_string;
      for (const auto& part : evt_(role)) {
        const unsigned short i = part.id();
        pyjets_.p[0][i] = part.momentum().px();
        pyjets_.p[1][i] = part.momentum().py();
        pyjets_.p[2][i] = part.momentum().pz();
        pyjets_.p[3][i] = part.momentum().energy();
        pyjets_.p[4][i] = part.momentum().mass();
        try {
          pyjets_.k[0][i] = pythia6Status((int)part.status());
        } catch (const std::out_of_range&) {
          evt_.dump();
          throw CG_FATAL("EventInterface") << "Failed to retrieve a Pythia 6 particle status translation for "
                                           << "CepGen status " << (int)part.status() << "!";
        }
        pyjets_.k[1][i] = part.integerPdgId();
        const auto& moth = part.mothers();
        pyjets_.k[2][i] = moth.empty() ? 0                   // no mother
                                       : *moth.begin() + 1;  // mother
        const auto& daug = part.daughters();
        if (daug.empty())  // no daughters
          pyjets_.k[3][i] = pyjets_.k[4][i] = 0;
        else {
          pyjets_.k[3][i] = *daug.begin() + 1;   // daughter 1
          pyjets_.k[4][i] = *daug.rbegin() + 1;  // daughter 2
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
      //--- at most one string per role
      if (!evt_string.empty())
        evt_strings_.emplace_back(evt_string);
    }

    //--- loop over the strings to bind everything together
    for (const auto& evt_string : evt_strings_) {
      if (evt_string.size() < 2)
        continue;

      CG_DEBUG_LOOP("EventInterface").log([&](auto& dbg) {
        dbg << "Joining " << utils::s("particle", evt_string.size()) << " with " << evt_[evt_string[0]].role()
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
    const auto old_npart = pyjets_.n;
    pyexec();
    //--- update the event
    for (int p = old_npart; p < pyjets_.n; ++p) {
      // filter the first particles already present in the event
      const auto pdg_id = std::abs(pyjets_.k[1][p]);
      checkPDGid(pdg_id);

      const auto moth_id = pyjets_.k[2][p] - 1;
      const auto role = pyjets_.k[2][p] != 0 ? evt_[moth_id].role()  // child particle inherits its mother's role
                                             : Particle::Role::UnknownRole;

      auto& pa = evt_.addParticle(role).get();
      pa.setId(p);
      pa.setStatus(cepgenStatus(pyjets_.k[0][p]));
      pa.setPdgId((long)pyjets_.k[1][p]);
      pa.setMomentum(
          Momentum(pyjets_.p[0][p], pyjets_.p[1][p], pyjets_.p[2][p], pyjets_.p[3][p]).setMass(pyjets_.p[4][p]));
      // define particle parentage
      auto& moth = evt_[moth_id];
      if (role != Particle::Role::UnknownRole)
        moth.setStatus(role == Particle::Role::CentralSystem ? Particle::Status::Resonance
                                                             : Particle::Status::Fragmented);
      pa.addMother(moth);
    }
  }

  std::pair<short, short> EventInterface::pickPartonsContent() {
    const auto ranudq = rnd_->uniformInt(0, 9);
    if (ranudq < 1)
      return {PDG::down, 2203};  // (d,uu1)
    if (ranudq < 5)
      return {PDG::up, 2101};  // (u,ud0)
    return {PDG::up, 2103};    // (u,ud1)
  }
}  // namespace cepgen::pythia6
