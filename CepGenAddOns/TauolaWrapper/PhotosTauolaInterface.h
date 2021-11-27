#ifndef CepGenAddOns_TauolaWrapper_PhotosTauolaInterface_h
#define CepGenAddOns_TauolaWrapper_PhotosTauolaInterface_h

#include "CepGen/Event/Event.h"
#include "CepGen/Physics/ParticleProperties.h"

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
      PhotosTauolaParticle(PhotosTauolaEvent<E, P>* event, const Particle& part);
      /// Remove all secondary particles instances created
      ~PhotosTauolaParticle();

      /// Create a new instance of a particle, disconnected from the event history
      PhotosTauolaParticle* createNewParticle(
          int pdg, int status, double mass, double px, double py, double pz, double e) override;
      /// Dump the particle attributes
      void print() override;

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
      void setMothers(std::vector<P*> mothers) override;
      /// Retrieve a list of parents from the event content
      std::vector<P*> getMothers() override;
      /// Specify a list of pointers to the secondary products
      void setDaughters(std::vector<P*> daughters) override;
      /// Retrieve a list of pointers to secondary products from the event content
      std::vector<P*> getDaughters() override;

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
      PhotosTauolaEvent(const Event& evt, const pdgid_t pdg = PDG::invalid);
      ~PhotosTauolaEvent();

      std::vector<P*> findParticles(int pdg) override;
      std::vector<P*> findStableParticles(int pdg) override;

    private:
      std::vector<P*> decay_particles_;
      pdgid_t spec_pdg_id_;
    };
  }  // namespace io
}  // namespace cepgen

#endif
