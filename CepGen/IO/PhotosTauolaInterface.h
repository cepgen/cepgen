#ifndef CepGen_IO_PhotosTauolaInterface_h
#define CepGen_IO_PhotosTauolaInterface_h

#include "CepGen/Core/Exception.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Event/Event.h"

namespace cepgen
{
  namespace io
  {
    template<typename E, typename P> class PhotosTauolaEvent;
    /// Interface to particles objects for Photos++ and Tauola++
    /// \tparam E Photos/Tauola event base object
    /// \tparam P Photos/Tauola particle base object
    template<typename E, typename P>
    class PhotosTauolaParticle : public P, public Particle
    {
      public:
        PhotosTauolaParticle() = default;
        inline PhotosTauolaParticle( PhotosTauolaEvent<E,P>* event, const Particle& part ) :
          Particle( part ), event_( event ) {
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
        inline ~PhotosTauolaParticle() {
          for ( auto& coll : { mothers_, daughters_, secondary_parts_ } )
            for ( auto& part : coll )
              delete part;
        }

        inline PhotosTauolaParticle* createNewParticle( int pdg, int status, double mass,
                                                        double px, double py, double pz, double e ) override {
          Particle part( Particle::Role::UnknownRole, pdg, (Particle::Status)status );
          part.setChargeSign( pdg/(unsigned int)pdg );
          part.setMomentum( Momentum::fromPxPyPzE( px, py, pz, e ) );
          part.setMass( mass );
          //event_->addParticle( Particle::CentralSystem ) = part;
          auto out = new PhotosTauolaParticle<E,P>( event_, part );
          secondary_parts_.emplace_back( out );
          return out;
        }
        inline void print() override { CG_INFO( "PhotosTauolaParticle" ) << *this; }

        void setBarcode( int id ) { id_ = id; }
        int getBarcode() override { return id_; }
        void setPdgID( int pdg ) override { setPdgId( (short)abs( pdg ) ); }
        int getPdgID() override { return integerPdgId(); }
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

        void setMothers( std::vector<P*> mothers ) override {
          for ( const auto& moth : mothers )
            addMother( *dynamic_cast<PhotosTauolaParticle*>( moth ) );
        }
        std::vector<P*> getMothers() override {
          if ( !mothers_.empty() )
            return mothers_;
          for ( const auto& moth : mothers() )
            mothers_.emplace_back( new PhotosTauolaParticle( event_, event_->operator[]( moth ) ) );
          return mothers_;
        }
        void setDaughters( std::vector<P*> daughters ) override {
          for ( const auto& daugh : daughters )
            addDaughter( *dynamic_cast<PhotosTauolaParticle*>( daugh ) );
        }
        std::vector<P*> getDaughters() override {
          if ( !daughters_.empty() )
            return daughters_;
          for ( const auto& daugh : daughters() )
            daughters_.emplace_back( new PhotosTauolaParticle( event_, event_->operator[]( daugh ) ) );
          return daughters_;
        }

      private:
        std::vector<P*> mothers_, daughters_;
        std::vector<P*> secondary_parts_;
        PhotosTauolaEvent<E,P>* event_;
    };

    /// Interface to events objects for Photos++ and Tauola++
    /// \tparam E Photos/Tauola event base object
    /// \tparam P Photos/Tauola particle base object
    template<typename E, typename P>
    class PhotosTauolaEvent : public E, public Event
    {
      public:
        inline PhotosTauolaEvent( const Event& evt, const pdgid_t pdg = PDG::invalid ) :
          Event( evt.compressed() ), spec_pdg_id_( pdg ) {
          //--- loop to add particles and associate parentage
          for ( size_t i = 0; i < size(); ++i ) {
            const auto& part = operator[]( i );
            particles_.emplace_back( new PhotosTauolaParticle<E,P>( this, part ) );
            const auto& moth = part.mothers(), &daugh = part.daughters();
            std::vector<P*> ml;
            for ( const auto& mt_id : moth )
              ml.emplace_back( particles_[mt_id] );
            particles_[i]->setMothers( ml );
          }
          //--- second loop to associate daughters (needs the full event to be filled)
          for ( size_t i = 0; i < size(); ++i ) {
            std::vector<P*> dl;
            for ( const auto& dg_id : operator[]( i ).daughters() )
              dl.emplace_back( particles_[dg_id] );
            particles_[i]->setDaughters( dl );
          }
        }
        inline ~PhotosTauolaEvent() {
          for ( size_t i = 0; i < particles_.size(); ++i )
            delete particles_[i];
        }

        inline std::vector<P*> findParticles( int pdg ) override {
          std::vector<P*> out;
          for ( auto& part : particles_ )
            if ( part->getPdgID() == pdg )
              out.emplace_back( part );
          return out;
        }
        inline std::vector<P*> findStableParticles( int pdg ) override {
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
