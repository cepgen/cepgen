#include "CepGen/Hadronisers/GenericHadroniser.h"
#include "CepGen/Hadronisers/HadronisersHandler.h"

#include "CepGen/IO/HepMCEventInterface.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"

#include <Tauola/Tauola.h>
#ifdef HEPMC3
#  include "CepGen/Event/Event.h"
#  include "CepGen/Core/utils.h"
#  include <Tauola/TauolaEvent.h>
#else
#  include <Tauola/TauolaHepMCEvent.h>
#endif
#include <Tauola/Log.h>

using namespace Tauolapp;

namespace cepgen
{
  namespace hadr
  {
    /// Interface to the Tauola decay routine
    class TauolaFilter : public GenericHadroniser
    {
      public:
        explicit TauolaFilter( const ParametersList& );
        ~TauolaFilter();

        void setParameters( const Parameters& ) override {}
        inline void readString( const char* param ) override {}
        void init() override;
        bool run( Event& ev, double& weight, bool full ) override;

        void setCrossSection( double xsec, double xsec_err ) override {}

      private:
        const ParametersList pol_states_, rad_states_;
#ifndef HEPMC3
        std::unique_ptr<HepMC::CepGenEvent> event_;
#endif
    };

#ifdef HEPMC3
    class CepGenTauolaEvent : public TauolaEvent, public Event
    {
      public:
        CepGenTauolaEvent( const Event& );
        ~CepGenTauolaEvent();

        std::vector<TauolaParticle*> findParticles( int ) override;
        std::vector<TauolaParticle*> findStableParticles( int ) override;

      private:
        std::vector<TauolaParticle*> particles_;
    };

    class CepGenTauolaParticle : public TauolaParticle, public Particle
    {
      public:
        CepGenTauolaParticle() = default;
        CepGenTauolaParticle( const Particle& );

        CepGenTauolaParticle* createNewParticle( int, int, double, double, double, double, double ) override;
        void print() override;

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

        void setMothers( std::vector<TauolaParticle*> moth ) override { mothers_ = moth; }
        std::vector<TauolaParticle*> getMothers() override { return mothers_; }
        void setDaughters( std::vector<TauolaParticle*> daugh ) override { daughters_ = daugh; }
        std::vector<TauolaParticle*> getDaughters() override { return daughters_; }

      private:
        std::vector<TauolaParticle*> mothers_, daughters_;
    };
#endif

    TauolaFilter::TauolaFilter( const ParametersList& params ) :
      GenericHadroniser( params, "tauola" ),
      pol_states_( params.get<ParametersList>( "polarisations" ) ),
      rad_states_( params.get<ParametersList>( "radiations" ) )
#ifndef HEPMC3
      , event_( new HepMC::CepGenEvent )
#endif
    {}

    TauolaFilter::~TauolaFilter()
    {
      Log::SummaryAtExit();
    }

    void
    TauolaFilter::init()
    {
      Tauola::setUnits( Tauola::GEV, Tauola::MM );
      //Tauola::setSeed( seed_ );
      Tauola::initialize();
      //--- spin correlations
      Tauola::spin_correlation.setAll( pol_states_.get<bool>( "full", true ) );
      Tauola::spin_correlation.GAMMA = pol_states_.get<bool>( "GAMMA", true );
      Tauola::spin_correlation.Z0 = pol_states_.get<bool>( "Z0", true );
      Tauola::spin_correlation.HIGGS = pol_states_.get<bool>( "HIGGS", true );
      Tauola::spin_correlation.HIGGS_H = pol_states_.get<bool>( "HIGGS_H", true );
      Tauola::spin_correlation.HIGGS_A = pol_states_.get<bool>( "HIGGS_A", true );
      Tauola::spin_correlation.HIGGS_PLUS = pol_states_.get<bool>( "HIGGS_PLUS", true );
      Tauola::spin_correlation.HIGGS_MINUS = pol_states_.get<bool>( "HIGGS_MINUS", true );
      Tauola::spin_correlation.W_PLUS = pol_states_.get<bool>( "W_PLUS", true );
      Tauola::spin_correlation.W_MINUS = pol_states_.get<bool>( "W_MINUS", true );
      //--- radiation states
      Tauola::setRadiation( rad_states_.get<bool>( "enable", true ) );
      const auto rad_cutoff = rad_states_.get<double>( "cutoff", -1. );
      if ( rad_cutoff > 0. )
        Tauola::setRadiationCutOff( rad_cutoff );
    }

    bool
    TauolaFilter::run( Event& ev, double& weight, bool full )
    {
      weight = 1.;

#ifndef HEPMC3
      event_->feedEvent( ev );
      event_->print();
      TauolaHepMCEvent evt( event_.get() );
#else
      CepGenTauolaEvent evt( ev );
#endif
      //evt.undecayTaus();
      evt.decayTaus();

      return true;
    }

#ifdef HEPMC3
    //----- Event interface

    CepGenTauolaEvent::CepGenTauolaEvent( const Event& evt ) :
      Event( evt )
    {
      //--- first loop to add particles
      for ( size_t i = 0; i < evt.size(); ++i )
        particles_.emplace_back( new CepGenTauolaParticle( evt[i] ) );
      //--- second loop to associate parentages
      for ( size_t i = 0; i < evt.size(); ++i ) {
        const auto& part = evt[i];
        const auto& daugh = part.daughters();
        if ( !daugh.empty() ) {
          std::vector<TauolaParticle*> dl;
          for ( const auto& dg_id : daugh )
            dl.emplace_back( particles_[dg_id] );
          particles_[i]->setDaughters( dl );
        }
        const auto& moth = part.mothers();
        if ( !moth.empty() ) {
          std::vector<TauolaParticle*> ml;
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

    CepGenTauolaEvent::~CepGenTauolaEvent()
    {
      for ( size_t i = 0; i < particles_.size(); ++i )
        delete particles_[i];
    }

    std::vector<TauolaParticle*>
    CepGenTauolaEvent::findParticles( int pdg )
    {
      std::vector<TauolaParticle*> out;
      for ( auto& part : particles_ )
        if ( part->getPdgID() == pdg )
          out.emplace_back( part );
      return out;
    }

    std::vector<TauolaParticle*>
    CepGenTauolaEvent::findStableParticles( int pdg )
    {
      std::vector<TauolaParticle*> out;
      for ( auto& part : findParticles( pdg ) )
        if ( part->getStatus() == TauolaParticle::STABLE )
          out.emplace_back( part );
      return out;
    }

    //----- Particle interface

    CepGenTauolaParticle::CepGenTauolaParticle( const Particle& part ) :
      Particle( part )
    {
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

    CepGenTauolaParticle*
    CepGenTauolaParticle::createNewParticle( int pdg, int status, double mass, double px, double py, double pz, double e )
    {
      Particle part( Particle::Role::UnknownRole, pdg, (Particle::Status)status );
      part.setChargeSign( pdg/(unsigned int)pdg );
      part.setMomentum( Particle::Momentum::fromPxPyPzE( px, py, pz, e ) );
      part.setMass( mass );
      return new CepGenTauolaParticle( part );
    }

    void
    CepGenTauolaParticle::print()
    {
      CG_INFO( "TauolaParticle" )
        << "TauolaParticle{pdg=" << getPdgID()
        << ",status=" << getStatus()
        << ",mom=(" << getPx() << "," << getPy() << "," << getPz() << ";" << getE() << ")"
        << ",m=" << getMass()
        << "," << utils::s( "mother", getMothers().size() )
        << "," << utils::s( "daughter", getDaughters().size() )
        << "}";
    }
#endif
  }
}

// register hadroniser and define alias
REGISTER_HADRONISER( tauola, TauolaFilter )

