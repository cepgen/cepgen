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
      auto out = new PhotosTauolaParticle<E, P>(event_, part);
      secondary_parts_.emplace_back(out);
      return out;
    }

    template <typename E, typename P>
    void PhotosTauolaParticle<E, P>::print() {
      CG_INFO("PhotosTauolaParticle:print") << *this;
    }

    template <typename E, typename P>
    void PhotosTauolaParticle<E, P>::setMothers(std::vector<P*> mothers) {
      for (const auto& moth : mothers) {
        auto&& part = dynamic_cast<PhotosTauolaParticle*>(moth);
        part->setStatus((int)Particle::Status::Propagator);
        Particle::addMother(*part);
      }
    }

    template <typename E, typename P>
    std::vector<P*> PhotosTauolaParticle<E, P>::getMothers() {
      if (!mothers_.empty())
        return mothers_;
      for (const auto& moth : mothers())
        if (moth >= 0)
          mothers_.emplace_back(new PhotosTauolaParticle(event_, event_->operator[](moth)));
      CG_DEBUG("PhotosTauolaParticle:getMothers").log([this](auto& log) {
        log << "List of mothers:";
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
        Particle::addDaughter(*part);
      }
    }

    template <typename E, typename P>
    std::vector<P*> PhotosTauolaParticle<E, P>::getDaughters() {
      if (daughters_.empty())
        for (const auto& daugh : daughters())
          if (daugh >= 0)
            daughters_.emplace_back(new PhotosTauolaParticle(event_, event_->operator[](daugh)));
      //CG_INFO("")<<getBarcode();for(const auto& p : daughters_)p->print();
      return daughters_;
    }

    //----- Event interface

    template <typename E, typename P>
    PhotosTauolaEvent<E, P>::PhotosTauolaEvent(const Event& evt, const pdgid_t pdg)
        : Event(evt.compress()), spec_pdg_id_(pdg) {}

    template <typename E, typename P>
    PhotosTauolaEvent<E, P>::~PhotosTauolaEvent() {
      for (size_t i = 0; i < decay_particles_.size(); ++i)
        delete decay_particles_[i];
    }

    template <typename E, typename P>
    std::vector<P*> PhotosTauolaEvent<E, P>::findParticles(int pdg) {
      //--- fill list of particles of interest if not already done
      if (decay_particles_.empty())
        for (auto& part : particles())
          if (abs(part.integerPdgId()) == pdg)
            decay_particles_.emplace_back(new PhotosTauolaParticle<E, P>(this, part));
      CG_DEBUG("PhotosTauolaEvent:findParticles").log([this](auto& log) {
        log << "Particles in event:";
        for (const auto* part : decay_particles_)
          log << "\n" << *dynamic_cast<const Particle*>(part);
      });
      return decay_particles_;
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
