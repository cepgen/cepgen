#include "CepGen/Hadronisers/Pythia8Hadroniser.h"
#include "CepGen/Hadronisers/PythiaEventInterface.h"
#include "CepGen/Hadronisers/HadronisersHandler.h"

#include "CepGen/Parameters.h"
#include "CepGen/Physics/Kinematics.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/PDG.h"

#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"

#include "CepGen/Event/Event.h"
#include "CepGen/Event/Particle.h"

#include "CepGen/Version.h"

namespace cepgen
{
  namespace hadr
  {
    Pythia8Hadroniser::Pythia8Hadroniser( const ParametersList& plist ) :
      GenericHadroniser( plist, "pythia8" ),
#ifdef PYTHIA8
      pythia_( new Pythia8::Pythia ), cg_evt_( new Pythia8::CepGenEvent ),
#endif
      full_evt_( false ), offset_( 0 ), first_evt_( true )
    {}

    void
    Pythia8Hadroniser::setParameters( const Parameters& params )
    {
      params_ = &params;
#ifdef PYTHIA8
      cg_evt_->initialise( params );
      pythia_->setLHAupPtr( (Pythia8::LHAup*)cg_evt_.get() );
      pythia_->settings.parm( "Beams:idA", (short)params_->kinematics.incoming_beams.first.pdg );
      pythia_->settings.parm( "Beams:idB", (short)params_->kinematics.incoming_beams.second.pdg );
      // specify we will be using a LHA input
      pythia_->settings.mode( "Beams:frameType", 5 );
      pythia_->settings.parm( "Beams:eCM", params_->kinematics.sqrtS() );
#endif
      for ( const auto& pdgid : params_->kinematics.minimum_final_state )
        min_ids_.emplace_back( (unsigned short)pdgid );
    }

    Pythia8Hadroniser::~Pythia8Hadroniser()
    {
#ifdef PYTHIA8
      pythia_->settings.writeFile( "last_pythia_config.cmd", false );
#endif
    }

    void
    Pythia8Hadroniser::readString( const char* param )
    {
#ifdef PYTHIA8
      if ( !pythia_->readString( param ) )
        throw CG_FATAL( "Pythia8Hadroniser" ) << "The Pythia8 core failed to parse the following setting:\n\t" << param;
#endif
    }

    void
    Pythia8Hadroniser::init()
    {
#ifdef PYTHIA8
      if ( pythia_->settings.flag( "ProcessLevel:all" ) != full_evt_ )
        pythia_->settings.flag( "ProcessLevel:all", full_evt_ );

      if ( seed_ == -1ll )
        pythia_->settings.flag( "Random:setSeed", false );
      else {
        pythia_->settings.flag( "Random:setSeed", true );
        pythia_->settings.mode( "Random:seed", seed_ );
      }

      switch ( params_->kinematics.mode ) {
        case KinematicsMode::ElasticElastic: {
          pythia_->settings.mode( "BeamRemnants:unresolvedHadron", 3 );
        } break;
        case KinematicsMode::InelasticElastic: {
          pythia_->settings.mode( "BeamRemnants:unresolvedHadron", 2 );
        } break;
        case KinematicsMode::ElasticInelastic: {
          pythia_->settings.mode( "BeamRemnants:unresolvedHadron", 1 );
        } break;
        case KinematicsMode::InelasticInelastic: default: {
          pythia_->settings.mode( "BeamRemnants:unresolvedHadron", 0 );
        } break;
      }

      if ( !pythia_->init() )
        throw CG_FATAL( "Pythia8Hadroniser" )
          << "Failed to initialise the Pythia8 core!\n\t"
          << "See the message above for more details.";
#else
      throw CG_FATAL( "Pythia8Hadroniser" )
        << "Pythia8 is not linked to this instance!";
#endif
    }

    void
    Pythia8Hadroniser::setCrossSection( double xsec, double xsec_err )
    {
#ifdef PYTHIA8
      cg_evt_->setCrossSection( 0, xsec, xsec_err );
#endif
    }

    bool
    Pythia8Hadroniser::run( Event& ev, double& weight, bool full )
    {
      //--- initialise the event weight before running any decay algorithm
      weight = 1.;

#ifdef PYTHIA8
      if ( !full && !pythia_->settings.flag( "ProcessLevel:resonanceDecays" ) )
        return true;

      //--- switch full <-> partial event
      if ( full != full_evt_ ) {
        full_evt_ = full;
        init();
      }

      //===========================================================================================
      // convert our event into a custom LHA format
      //===========================================================================================

      cg_evt_->feedEvent( ev, full );
      //if ( full ) cg_evt_->listEvent();

      //===========================================================================================
      // launch the hadronisation / resonances decays, and update the event accordingly
      //===========================================================================================

      ev.num_hadronisation_trials = 0;
      while ( true ) {
        if ( ev.num_hadronisation_trials++ > max_trials_ )
          return false;
        //--- run the hadronisation/fragmentation algorithm
        if ( pythia_->next() ) {
          //--- hadronisation successful
          if ( first_evt_ && full ) {
            offset_ = 0;
            for ( unsigned short i = 1; i < pythia_->event.size(); ++i )
              if ( pythia_->event[i].status() == -12 ) // skip the incoming particles
                offset_++;
            first_evt_ = false;
          }
          break;
        }
      }

      //===========================================================================================
      // update the event content with Pythia's output
      //===========================================================================================

      updateEvent( ev, weight );
      return true;
#else
      throw CG_FATAL( "Pythia8Hadroniser" ) << "Pythia8 is not linked to this instance!";
#endif
    }

#ifdef PYTHIA8
    Particle&
    Pythia8Hadroniser::addParticle( Event& ev, const Pythia8::Particle& py_part, const Pythia8::Vec4& mom, unsigned short role ) const
    {
      Particle& op = ev.addParticle( (Particle::Role)role );
      op.setPdgId( static_cast<PDG>( abs( py_part.id() ) ), (short)( py_part.charge()/fabs( py_part.charge() ) ) );
      op.setStatus( py_part.isFinal()
        ? Particle::Status::FinalState
        : Particle::Status::Propagator );
      op.setMomentum( Particle::Momentum( mom.px(), mom.py(), mom.pz(), mom.e() ) );
      op.setMass( mom.mCalc() );
      cg_evt_->addCorresp( py_part.index()-offset_, op.id() );
      return op;
    }

    void
    Pythia8Hadroniser::updateEvent( Event& ev, double& weight ) const
    {
      for ( unsigned short i = 1+offset_; i < pythia_->event.size(); ++i ) {
        const Pythia8::Particle& p = pythia_->event[i];
        const unsigned short cg_id = cg_evt_->cepgenId( i-offset_ );
        if ( cg_id != Pythia8::CepGenEvent::INVALID_ID ) {
          //----- particle already in the event
          Particle& cg_part = ev[cg_id];
          //--- fragmentation result
          if ( cg_part.role() == Particle::OutgoingBeam1
            || cg_part.role() == Particle::OutgoingBeam2 ) {
            cg_part.setStatus( Particle::Status::Fragmented );
            continue;
          }
          //--- particle is not what we expect
          if ( p.idAbs() != abs( cg_part.integerPdgId() ) ) {
            CG_INFO( "Pythia8Hadroniser:update" ) << "LHAEVT event content:";
            cg_evt_->listEvent();
            CG_INFO( "Pythia8Hadroniser:update" ) << "Pythia event content:";
            pythia_->event.list();
            CG_INFO( "Pythia8Hadroniser:update" ) << "CepGen event content:";
            ev.dump();
            CG_INFO( "Pythia8Hadroniser:update" ) << "Correspondence:";
            cg_evt_->dumpCorresp();

            throw CG_FATAL( "Pythia8Hadroniser:update" )
              << "Event list corruption detected for (Pythia/CepGen) particle " << i << "/" << cg_id << ":\n\t"
              << "should be " << abs( p.id() ) << ", "
              << "got " << cg_part.integerPdgId() << "!";
          }
          //--- resonance decayed; apply branching ratio for this decay
          if ( p.status() < 0 && pythia_->settings.flag( "ProcessLevel:resonanceDecays" ) ) {
            for ( unsigned short j = 0; j < p.particleDataEntry().sizeChannels(); ++j )
              for ( const auto& d : p.daughterListRecursive() )
                if ( p.particleDataEntry().channel( j ).contains( pythia_->event[d].id() ) ) {
                  weight *= p.particleDataEntry().channel( j ).bRatio();
                  break;
                }
            //weight *= p.particleDataEntry().pickChannel().bRatio();
            cg_part.setStatus( Particle::Status::Resonance );
          }
        }
        else {
          //----- new particle to be added
          const unsigned short role = findRole( ev, p );
          switch ( (Particle::Role)role ) {
            default: break;
            case Particle::OutgoingBeam1: {
              ev[Particle::OutgoingBeam1][0].setStatus( Particle::Status::Fragmented );
              if ( abs( p.status() ) != 61 )
                break;
            } // no break!
            case Particle::OutgoingBeam2: {
              ev[Particle::OutgoingBeam2][0].setStatus( Particle::Status::Fragmented );
              if ( abs( p.status() ) != 61 )
                break;
            } // no break!
          }
          // found the role ; now we can add the particle
          Particle& cg_part = addParticle( ev, p, p.p(), role );
          for ( const auto& moth_id : p.motherList() ) {
            if ( moth_id <= offset_ )
              continue;
            const unsigned short moth_cg_id = cg_evt_->cepgenId( moth_id-offset_ );
            if ( moth_cg_id != Pythia8::CepGenEvent::INVALID_ID )
              cg_part.addMother( ev[moth_cg_id] );
            else
              cg_part.addMother( addParticle( ev, pythia_->event[moth_id], p.p(), role ) );
            if ( !p.isFinal() ) {
              if ( p.isResonance() || !p.daughterList().empty() )
                cg_part.setStatus( Particle::Status::Resonance );
              else
                cg_part.setStatus( Particle::Status::Undefined );
            }
          }
        }
      }
    }

    unsigned short
    Pythia8Hadroniser::findRole( const Event& ev, const Pythia8::Particle& p ) const
    {
      for ( const auto& par_id : p.motherList() ) {
        if ( par_id == 1 && offset_ > 0 )
          return (unsigned short)Particle::OutgoingBeam1;
        if ( par_id == 2 && offset_ > 0 )
          return (unsigned short)Particle::OutgoingBeam2;
        const unsigned short par_cg_id = cg_evt_->cepgenId( par_id-offset_ );
        if ( par_cg_id != Pythia8::CepGenEvent::INVALID_ID )
          return (unsigned short)ev[par_cg_id].role();
        return findRole( ev, pythia_->event[par_id] );
      }
      return (unsigned short)Particle::UnknownRole;
    }
#endif
  }
}
// register hadroniser and define alias
REGISTER_HADRONISER( pythia8, Pythia8Hadroniser )
