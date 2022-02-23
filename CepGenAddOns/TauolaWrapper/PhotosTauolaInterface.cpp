/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2021-2022  Laurent Forthomme
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

#include <Tauola/TauolaEvent.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/String.h"
#include "CepGenAddOns/TauolaWrapper/PhotosTauolaInterface.h"

namespace cepgen {
  namespace io {

    //----- Particle interface

    template <typename E, typename P>
    PhotosTauolaParticle<E, P>::PhotosTauolaParticle(PhotosTauolaEvent<E, P>* event, const Particle& part)
        : Particle(part), event_(event) {}

    template <typename E, typename P>
    PhotosTauolaParticle<E, P>::~PhotosTauolaParticle() {
      //--- clean all collections
      for (auto& coll : {mothers_, daughters_, secondary_parts_})
        for (size_t i = 0; i < coll.size(); ++i)
          delete coll[i];
    }

    template <typename E, typename P>
    PhotosTauolaParticle<E, P>* PhotosTauolaParticle<E, P>::createNewParticle(
        int pdg, int status, double mass, double px, double py, double pz, double e) {
      Particle part(Particle::Role::CentralSystem, pdg, (Particle::Status)status);
      part.setChargeSign(pdg / abs(pdg));
      part.setMomentum(Momentum::fromPxPyPzE(px, py, pz, e));
      part.setMass(mass);
      secondary_parts_.emplace_back(new PhotosTauolaParticle<E, P>(event_, part));
      CG_DEBUG("PhotosTauolaParticle:createNewParticle") << "New particle built: " << part << ".";
      return dynamic_cast<PhotosTauolaParticle<E, P>*>(*secondary_parts_.rbegin());
    }

    template <typename E, typename P>
    int PhotosTauolaParticle<E, P>::getPdgID() {
      CG_WARNING("") << integerPdgId();
      return Particle::integerPdgId();
    }

    template <typename E, typename P>
    void PhotosTauolaParticle<E, P>::setStatus(int status) {
      status_ = status;
    }

    template <typename E, typename P>
    void PhotosTauolaParticle<E, P>::print() {
      CG_INFO("PhotosTauolaParticle:print") << *this;
    }

    template <typename E, typename P>
    void PhotosTauolaParticle<E, P>::undecay() {
      CG_WARNING("");
    }

    template <typename E, typename P>
    void PhotosTauolaParticle<E, P>::checkMomentumConservation() {
      CG_WARNING("");
    }

    template <typename E, typename P>
    void PhotosTauolaParticle<E, P>::decayEndgame() {
      CG_WARNING("");
    }

    template <typename E, typename P>
    void PhotosTauolaParticle<E, P>::setMothers(std::vector<P*> mothers) {
      CG_WARNING("");
      for (const auto& moth : mothers) {
        auto&& part = dynamic_cast<PhotosTauolaParticle*>(moth);
        part->setStatus((int)Particle::Status::Propagator);
        mothers_.emplace_back(moth);
        Particle::addMother(*part);
      }
      CG_DEBUG("PhotosTauolaParticle:setMothers") << "New list of mothers: " << mothers_ << ".";
    }

    template <typename E, typename P>
    std::vector<P*> PhotosTauolaParticle<E, P>::getMothers() {
      if (!mothers_.empty())
        return mothers_;
      for (const auto& moth : mothers())
        if (moth >= 0)
          mothers_.emplace_back(new PhotosTauolaParticle(event_, event_->operator[](moth)));
      CG_DEBUG("PhotosTauolaParticle:getMothers").log([this](auto& log) {
        if (mothers_.empty()) {
          log << "No mothers for " << *this << ".";
          return;
        }
        log << "Mothers for " << *this << ":";
        for (const auto* moth : mothers_)
          log << "\n" << *dynamic_cast<const Particle*>(moth);
      });
      return mothers_;
    }

    template <typename E, typename P>
    void PhotosTauolaParticle<E, P>::setDaughters(std::vector<P*> daughters) {
      for (const auto& daugh : daughters) {
        auto&& part = dynamic_cast<PhotosTauolaParticle*>(daugh);
        part->setRole(role());  // child inherits its mother's role
        daughters_.emplace_back(daugh);
        Particle::addDaughter(*part);
      }
      CG_DEBUG("PhotosTauolaParticle:setDaughters") << "New list of daughters:" << daughters_ << ".";
    }

    template <typename E, typename P>
    std::vector<P*> PhotosTauolaParticle<E, P>::getDaughters() {
      if (daughters_.empty())
        for (const auto& daugh : daughters())
          if (daugh >= 0)
            daughters_.emplace_back(new PhotosTauolaParticle(event_, event_->operator[](daugh)));
      CG_DEBUG("PhotosTauolaParticle:getDaughters").log([this](auto& log) {
        if (daughters_.empty()) {
          log << "No daughters for " << *this << ".";
          return;
        }
        log << "Daughters for " << *this << ":";
        for (const auto* daugh : daughters_)
          log << "\n" << *dynamic_cast<const Particle*>(daugh);
      });
      return daughters_;
    }

    //----- Event interface

    template <typename E, typename P>
    PhotosTauolaEvent<E, P>::PhotosTauolaEvent(const Event& evt, const pdgid_t pdg)
        : Event(evt.compress()), spec_pdg_id_(pdg) {}

    template <typename E, typename P>
    PhotosTauolaEvent<E, P>::~PhotosTauolaEvent() {}

    template <typename E, typename P>
    void PhotosTauolaEvent<E, P>::eventEndgame() {
      CG_WARNING("");
    }

    template <typename E, typename P>
    std::vector<P*> PhotosTauolaEvent<E, P>::findParticles(int pdg) {
      std::vector<P*> out;
      //--- fill list of particles of interest if not already done
      for (auto& part : particles())
        if (abs(part.integerPdgId()) == pdg)
          out.emplace_back(new PhotosTauolaParticle<E, P>(this, part));
      CG_DEBUG("PhotosTauolaEvent:findParticles").log([&out](auto& log) {
        log << "Particles in event:";
        for (const auto* part : out)
          log << "\n" << *dynamic_cast<const Particle*>(part);
      });
      return out;
    }

    template <typename E, typename P>
    std::vector<P*> PhotosTauolaEvent<E, P>::findStableParticles(int pdg) {
      std::vector<P*> out;
      for (auto& part : findParticles(pdg)) {
        if (!part->hasDaughters())
          out.emplace_back(part);
        else {
          const auto& daugh = part->getDaughters();
          if (daugh.size() == 1)
            continue;  // weird parentage, particle will not be decayed
          if (daugh.size() == 2 &&
              (abs(daugh.at(0)->getPdgID()) == spec_pdg_id_ || abs(daugh.at(1)->getPdgID()) == spec_pdg_id_))
            continue;  // already decayed into a pair of particles of interest; skip it
          CG_WARNING("PhotosTauolaEvent") << "Particle with pdg code " << part->getPdgID() << " has already "
                                          << utils::s("daughter", daugh.size(), true) << ".";
        }
      }
      CG_DEBUG("PhotosTauolaEvent:findStableParticles").log([&out](auto& log) {
        log << "Stable particles in event:";
        for (const auto* part : out)
          log << "\n" << *dynamic_cast<const Particle*>(part);
      });
      return out;
    }

    template class PhotosTauolaEvent<Tauolapp::TauolaEvent, Tauolapp::TauolaParticle>;
  }  // namespace io
}  // namespace cepgen
