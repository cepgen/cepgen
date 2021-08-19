#ifndef CepGenAddOns_TauolaWrapper_PhotosTauolaInterface_h
#define CepGenAddOns_TauolaWrapper_PhotosTauolaInterface_h

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Physics/PDG.h"

namespace cepgen {
  namespace io {
    // Forward declaration
    template <typename E, typename P>
    class PhotosTauolaEvent;
    /// Interface to particles objects for Photos++ and Tauola++
    /// \tparam E Photos/Tauola event base object
    /// \tparam P Photos/Tauola particle base object
    template <typename E, typename P>
    class PhotosTauolaParticle : public P, public Particle {
    public:
      PhotosTauolaParticle() = default;
      inline PhotosTauolaParticle(PhotosTauolaEvent<E, P>* event, const Particle& part)
          : Particle(part), event_(event) {}
      /// Remove all secondary particles instances created
      inline ~PhotosTauolaParticle() {
        //--- clean all collections
        for (auto& coll : {mothers_, daughters_, secondary_parts_})
          for (size_t i = 0; i < coll.size(); ++i)
            delete coll[i];
      }

      /// Create a new instance of a particle, disconnected from the event history
      inline PhotosTauolaParticle* createNewParticle(
          int pdg, int status, double mass, double px, double py, double pz, double e) override {
        Particle part(Particle::Role::CentralSystem, pdg, (Particle::Status)status);
        part.setChargeSign(pdg / abs(pdg));
        part.setMomentum(Momentum::fromPxPyPzE(px, py, pz, e));
        part.setMass(mass);
        auto out = new PhotosTauolaParticle<E, P>(event_, part);
        secondary_parts_.emplace_back(out);
        return out;
      }
      /// Dump the particle attributes
      inline void print() override { throw CG_FATAL("PhotosTauolaParticle") << *this; }

      /// Specify the particle unique identifier
      void setBarcode(int id) { id_ = id; }
      /// Particle unique identifier in the event
      int getBarcode() override { return id_; }
      /// Set the particle ID
      void setPdgID(int pdg) override { Particle::setPdgId((long)pdg); }
      /// Particle ID
      int getPdgID() override { return Particle::integerPdgId(); }
      void setStatus(int status) override { status_ = status; }
      /// Particle status
      int getStatus() override { return status_; }
      void setPx(double px) override { momentum_.setPx(px); }
      /// Horizontal component of the momentum
      double getPx() override { return momentum_.px(); }
      void setPy(double py) override { momentum_.setPy(py); }
      /// Vertical component of the momentum
      double getPy() override { return momentum_.py(); }
      void setPz(double pz) override { momentum_.setPz(pz); }
      /// Longitudinal component of the momentum
      double getPz() override { return momentum_.pz(); }
      void setE(double e) override { setE(e); }
      /// Particle energy
      double getE() override { return momentum_.energy(); }
      void setMass(double m) override { mass_ = m; }
      /// Specify a list of pointers to the parents
      void setMothers(std::vector<P*> mothers) override {
        for (const auto& moth : mothers) {
          auto&& part = dynamic_cast<PhotosTauolaParticle*>(moth);
          part->setStatus((int)Particle::Status::Propagator);
          Particle::addMother(*part);
        }
      }
      /// Retrieve a list of parents from the event content
      std::vector<P*> getMothers() override {
        if (!mothers_.empty())
          return mothers_;
        for (const auto& moth : mothers())
          if (moth >= 0)
            mothers_.emplace_back(new PhotosTauolaParticle(event_, event_->operator[](moth)));
        return mothers_;
      }
      /// Specify a list of pointers to the secondary products
      void setDaughters(std::vector<P*> daughters) override {
        for (const auto& daugh : daughters) {
          auto&& part = dynamic_cast<PhotosTauolaParticle*>(daugh);
          part->setRole(role());  // child inherits its mother's role
          Particle::addDaughter(*part);
        }
      }
      /// Retrieve a list of pointers to secondary products from the event content
      std::vector<P*> getDaughters() override {
        if (daughters_.empty())
          for (const auto& daugh : daughters())
            if (daugh >= 0)
              daughters_.emplace_back(new PhotosTauolaParticle(event_, event_->operator[](daugh)));
        //CG_INFO("")<<getBarcode();for(const auto& p : daughters_)p->print();
        return daughters_;
      }

    private:
      std::vector<P*> mothers_, daughters_;
      std::vector<P*> secondary_parts_;
      PhotosTauolaEvent<E, P>* event_;  // non-owning, only treated as reference
    };

    /// Interface to events objects for Photos++ and Tauola++
    /// \tparam E Photos/Tauola event base object
    /// \tparam P Photos/Tauola particle base object
    template <typename E, typename P>
    class PhotosTauolaEvent : public E, public Event {
    public:
      inline PhotosTauolaEvent(const Event& evt, const pdgid_t pdg = PDG::invalid)
          : Event(evt.compress()), spec_pdg_id_(pdg) {}
      inline ~PhotosTauolaEvent() {
        for (size_t i = 0; i < decay_particles_.size(); ++i)
          delete decay_particles_[i];
      }

      inline std::vector<P*> findParticles(int pdg) override {
        //--- fill list of particles of interest if not already done
        if (decay_particles_.empty())
          for (auto& part : particles())
            if (abs(part.integerPdgId()) == pdg)
              decay_particles_.emplace_back(new PhotosTauolaParticle<E, P>(this, part));
        return decay_particles_;
      }
      inline std::vector<P*> findStableParticles(int pdg) override {
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
        return out;
      }

    private:
      std::vector<P*> decay_particles_;
      pdgid_t spec_pdg_id_;
    };
  }  // namespace io
}  // namespace cepgen

#endif
