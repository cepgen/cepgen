/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2023  Laurent Forthomme
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
#include "CepGen/Utils/String.h"
#include "CepGenAddOns/Pythia6Wrapper/EventInterface.h"
#include "CepGenAddOns/Pythia6Wrapper/Pythia6Interface.h"

namespace pythia6 {
  EventInterface::EventInterface(cepgen::Event& event, const cepgen::mode::Kinematics kin_mode)
      : evt_(event), rnd_phi_(0., 2. * M_PI), rnd_cos_theta_(-1., 1.), rnd_qdq_(0., 9.) {
    if (kin_mode == cepgen::mode::Kinematics::InelasticElastic ||
        kin_mode == cepgen::mode::Kinematics::InelasticInelastic)
      roles_.emplace_back(cepgen::Particle::Role::OutgoingBeam1);
    if (kin_mode == cepgen::mode::Kinematics::ElasticInelastic ||
        kin_mode == cepgen::mode::Kinematics::InelasticInelastic)
      roles_.emplace_back(cepgen::Particle::Role::OutgoingBeam2);
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
      const double phi = rnd_phi_(rnd_gen_), theta = acos(rnd_cos_theta_(rnd_gen_));

      // compute momentum of decay particles from MX
      const double px2 = std::pow(cepgen::utils::energyFromW(part.momentum().mass(), mdq2, mq2), 2) - mq2;
      if (px2 < 0.) {
        CG_WARNING("EventInterface:prepareHadronisation") << "Invalid remnants kinematics for " << role << ".";
        continue;
      }
      const double px = std::sqrt(px2);

      //--- build 4-vectors and boost decay particles
      auto pdq = cepgen::Momentum::fromPThetaPhiE(px, theta, phi, std::hypot(px, mdq));
      auto pq = -pdq;
      pq.setEnergy(std::hypot(px, mq));

      const auto part_id = part.id();

      // singlet
      auto& quark = evt_.addParticle(role).get();
      quark.addMother(evt_[part_id]);
      quark.setPdgId(partons.first, +1);
      quark.setStatus(cepgen::Particle::Status::Unfragmented);
      quark.setMomentum(pq.lorentzBoost(evt_[part_id].momentum()));

      // quark doublet
      auto& diquark = evt_.addParticle(role).get();
      diquark.addMother(evt_[part_id]);
      diquark.setPdgId(partons.second, +1);
      diquark.setStatus(cepgen::Particle::Status::Unfragmented);
      diquark.setMomentum(pdq.lorentzBoost(evt_[part_id].momentum()));

      evt_[part_id].setStatus(cepgen::Particle::Status::Fragmented);
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

        if (part.status() == cepgen::Particle::Status::Unfragmented) {
          pyjets_.k[0][i] = 1;  // PYTHIA/JETSET workaround
          evt_string.emplace_back(part.id() + 1);
        } else if (part.status() == cepgen::Particle::Status::Undecayed)
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
        dbg << "Joining " << cepgen::utils::s("particle", evt_string.size()) << " with " << evt_[evt_string[0]].role()
            << " role"
            << " in a same string";
        for (const auto& part_id : evt_string) {
          if (part_id != -1)
            dbg << cepgen::utils::format("\n\t * %2d (pdgId=%4d)", part_id, pyjets_.k[1][part_id - 1]);
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
                                             : cepgen::Particle::Role::UnknownRole;

      auto& pa = evt_.addParticle(role).get();
      pa.setId(p);
      pa.setStatus(cepgenStatus(pyjets_.k[0][p]));
      pa.setPdgId((long)pyjets_.k[1][p]);
      pa.setMomentum(
          cepgen::Momentum(pyjets_.p[0][p], pyjets_.p[1][p], pyjets_.p[2][p], pyjets_.p[3][p]).setMass(pyjets_.p[4][p]));
      // define particle parentage
      auto& moth = evt_[moth_id];
      if (role != cepgen::Particle::Role::UnknownRole)
        moth.setStatus(role == cepgen::Particle::Role::CentralSystem ? cepgen::Particle::Status::Resonance
                                                                     : cepgen::Particle::Status::Fragmented);
      pa.addMother(moth);
    }
  }

  std::pair<short, short> EventInterface::pickPartonsContent() {
    const auto ranudq = rnd_qdq_(rnd_gen_);
    if (ranudq < 1.)
      return {cepgen::PDG::down, 2203};  // (d,uu1)
    if (ranudq < 5.)
      return {cepgen::PDG::up, 2101};  // (u,ud0)
    return {cepgen::PDG::up, 2103};    // (u,ud1)
  }
}  // namespace pythia6
