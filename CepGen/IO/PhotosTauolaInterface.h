#ifndef CepGen_IO_PhotosTauolaInterface_h
#define CepGen_IO_PhotosTauolaInterface_h

#include "CepGen/Core/Exception.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Event/Event.h"

namespace cepgen
{
  namespace io
  {
    template<typename P>
    class PhotosTauolaParticle : public P, public Particle
    {
      public:
        PhotosTauolaParticle() = default;
        PhotosTauolaParticle( const Particle& part ) :
          Particle( part ) {
          setBarcode( part.id() );
          setPdgID( part.integerPdgId() );
          setStatus( (int)part.status() );
          const auto& mom = part.momentum();
          setPx( mom.px() );
          setPy( mom.py() );
          setPz( mom.pz() );
          setEnergy( mom.energy() );
          setMass( part.mass() );
        }

        PhotosTauolaParticle* createNewParticle( int pdg, int status, double mass, double px, double py, double pz, double e ) override {
          Particle part( Particle::Role::UnknownRole, pdg, (Particle::Status)status );
          part.setChargeSign( pdg/(unsigned int)pdg );
          part.setMomentum( Momentum::fromPxPyPzE( px, py, pz, e ) );
          part.setMass( mass );
          return new PhotosTauolaParticle<P>( part );
        }
        void print() override { CG_INFO( "PhotosTauolaParticle" ) << *this; }

        void setBarcode( int id ) { id_ = id; }
        int getBarcode() override { return id_; }
        void setPdgID( int pdg ) override { pdg_id_ = (pdgid_t)pdg; }
        int getPdgID() override { return (int)pdg_id_; }
        void setStatus( int status ) override { status_ = (Particle::Status)status; }
        int getStatus() override { return (int)status_; }

        void setPx( double px ) override { momentum_[0] = px; }
        double getPx() override { return momentum_.px(); }
        void setPy( double py ) override { momentum_[1] = py; }
        double getPy() override { return momentum_.py(); }
        void setPz( double pz ) override { momentum_[2] = pz; }
        double getPz() override { return momentum_.pz(); }
        void setE( double e ) override { momentum_[3] = e; }
        double getE() override { return energy(); }
        void setMass( double m ) override { mass_ = m; }

        void setMothers( std::vector<P*> moth ) override { mothers_ = moth; }
        std::vector<P*> getMothers() override { return mothers_; }
        void setDaughters( std::vector<P*> daugh ) override { daughters_ = daugh; }
        std::vector<P*> getDaughters() override { return daughters_; }

      private:
        std::vector<P*> mothers_, daughters_;
    };

    template<typename E, typename P>
    class PhotosTauolaEvent : public E, public Event
    {
      public:
        inline PhotosTauolaEvent( const Event& evt, const pdgid_t pdg = PDG::invalid ) :
          Event( evt ), spec_pdg_id_( pdg ) {
          //--- first loop to add particles
          for ( size_t i = 0; i < evt.size(); ++i )
            particles_.emplace_back( new PhotosTauolaParticle<P>( evt[i] ) );
          //--- second loop to associate parentages
          for ( size_t i = 0; i < evt.size(); ++i ) {
            const auto& part = evt[i];
            const auto& daugh = part.daughters();
            if ( !daugh.empty() ) {
              std::vector<P*> dl;
              for ( const auto& dg_id : daugh )
                dl.emplace_back( particles_[dg_id] );
              particles_[i]->setDaughters( dl );
            }
            const auto& moth = part.mothers();
            if ( !moth.empty() ) {
              std::vector<P*> ml;
              for ( const auto& mt_id : moth ) {
                ml.emplace_back( particles_[mt_id] );
                if ( part.role() == Particle::Role::Intermediate )
                  particles_[i]->setPdgID( particles_[mt_id]->getPdgID() );
              }
              particles_[i]->setMothers( ml );
            }
            particles_[i]->print();
          }
        }
        inline ~PhotosTauolaEvent() {
          for ( size_t i = 0; i < particles_.size(); ++i )
            delete particles_[i];
        }

        std::vector<P*> findParticles( int pdg ) override {
          std::vector<P*> out;
          for ( auto& part : particles_ )
            if ( part->getPdgID() == pdg )
              out.emplace_back( part );
          return out;
        }
        std::vector<P*> findStableParticles( int pdg ) override {
          std::vector<P*> out;
          for ( auto& part : findParticles( pdg ) ) {
            if ( !part->hasDaughters() )
              out.emplace_back( part );
            else {
              const auto& daugh = part->getDaughters();
              if ( daugh.size() == 1 )
                continue;
              if ( daugh.size() == 2 && ( abs( daugh.at( 0 )->getPdgID() ) == spec_pdg_id_
                                       || abs( daugh.at( 1 )->getPdgID() ) == spec_pdg_id_ ) )
                continue;
              CG_WARNING( "PhotosTauolaEvent" ) << "Particle with pdg code " << part->getPdgID()
                <<" has already daughters.";
            }
          }
          return out;
        }

      private:
        std::vector<P*> particles_;
        pdgid_t spec_pdg_id_;
    };
  }
}

#endif
