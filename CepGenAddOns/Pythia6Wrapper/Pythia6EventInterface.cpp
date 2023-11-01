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
#include "CepGen/Utils/String.h"
#include "CepGenAddOns/Pythia6Wrapper/Pythia6EventInterface.h"
#include "CepGenAddOns/Pythia6Wrapper/Pythia6Interface.h"

namespace cepgen {
  namespace hadr {
    Pythia6EventInterface::Pythia6EventInterface()
        : rnd_phi_(0., 2. * M_PI), rnd_cos_theta_(-1., 1.), rnd_qdq_(0., 9.) {}

    void Pythia6EventInterface::feedEvent(const cepgen::Event& evt) {
      evt_ = evt;
      //--- initialising the string fragmentation variables
      evt_strings_.clear();

      pyjets_.n = 0;  // reinitialise the event content

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
            pyjets_.k[0][i] = pythia6::status((int)part.status());
          } catch (const std::out_of_range&) {
            evt_.dump();
            throw CG_FATAL("Pythia6Hadroniser") << "Failed to retrieve a Pythia 6 particle status translation for "
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

        CG_DEBUG_LOOP("Pythia6Hadroniser").log([&](auto& dbg) {
          dbg << "Joining " << utils::s("particle", evt_string.size()) << " with " << evt_[evt_string[0]].role()
              << " role"
              << " in a same string";
          for (const auto& part_id : evt_string) {
            if (part_id != -1)
              dbg << utils::format("\n\t * %2d (pdgId=%4d)", part_id, pyjets_.k[1][part_id - 1]);
          }
        });
        pythia6::pyjoin(evt_string);
      }

      CG_DEBUG_LOOP("Pythia6Hadroniser") << "Hadronisation preparation called.";

      for (const auto& part : evt_.particles()) {
        if (part.status() != Particle::Status::Unfragmented)
          continue;
        //--- only loop over all protons to be fragmented

        const auto partons = pickPartonsContent();
        const double mx2 = part.momentum().mass2();
        const double mq = pythia6::pymass(partons.first), mq2 = mq * mq;
        const double mdq = pythia6::pymass(partons.second), mdq2 = mdq * mdq;

        //--- choose random direction in MX frame
        const double phi = rnd_phi_(rnd_gen_), theta = acos(rnd_cos_theta_(rnd_gen_));  // theta angle

        //--- compute momentum of decay particles from MX
        const double px2 = 0.25 * std::pow(mx2 - mdq2 + mq2, 2) / mx2 - mq2;
        if (px2 < 0.) {
          CG_WARNING("Pythia6Hadroniser") << "Invalid remnants kinematics for " << part.role() << ".";
          return;
        }
        const double px = std::sqrt(px2);

        //--- build 4-vectors and boost decay particles
        auto pdq = Momentum::fromPThetaPhiE(px, theta, phi, std::hypot(px, mdq));
        auto pq = -pdq;
        pq.setEnergy(std::hypot(px, mq));

        //--- singlet
        auto& quark = evt_.addParticle(part.role()).get();
        quark.addMother(evt_[part.id()]);
        quark.setPdgId(partons.first, +1);
        quark.setStatus(Particle::Status::Unfragmented);
        quark.setMomentum(pq.lorentzBoost(part.momentum()));

        //--- doublet
        auto& diquark = evt_.addParticle(part.role()).get();
        diquark.addMother(evt_[part.id()]);
        diquark.setPdgId(partons.second, +1);
        diquark.setStatus(Particle::Status::Unfragmented);
        diquark.setMomentum(pdq.lorentzBoost(part.momentum()));

        evt_[part.id()].setStatus(Particle::Status::Fragmented);
      }
    }

    void Pythia6EventInterface::run() const { pythia6::pyexec(); }

    std::pair<short, short> Pythia6EventInterface::pickPartonsContent() const {
      const auto ranudq = rnd_qdq_(rnd_gen_);
      if (ranudq < 1.)
        return {PDG::down, 2203};  // (d,uu1)
      if (ranudq < 5.)
        return {PDG::up, 2101};  // (u,ud0)
      return {PDG::up, 2103};    // (u,ud1)
    }
  }  // namespace hadr
}  // namespace cepgen
