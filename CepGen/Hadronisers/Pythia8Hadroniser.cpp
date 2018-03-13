#ifdef PYTHIA8
#include "Pythia8Hadroniser.h"

#include "CepGen/Parameters.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Event/Particle.h"
#include "CepGen/Version.h"

namespace CepGen
{
  namespace Hadroniser
  {
    Pythia8Hadroniser::Pythia8Hadroniser( const Parameters& params ) :
      GenericHadroniser( "pythia8" ),
      max_attempts_( params.hadroniser_max_trials ),
      lhaevt_( new LHAEvent )
    {
#ifdef PYTHIA8
      pythia_.reset( new Pythia8::Pythia );
      pythia_->settings.parm( "Beams:idA", params.kinematics.inpdg.first );
      pythia_->settings.parm( "Beams:idB", params.kinematics.inpdg.second );
      pythia_->settings.parm( "Beams:eCM", params.kinematics.sqrtS() );
#endif
      for ( const auto& pdgid : params.kinematics.minimum_final_state )
        min_ids_.emplace_back( (unsigned short)pdgid );
    }

    Pythia8Hadroniser::~Pythia8Hadroniser()
    {}

    bool
    Pythia8Hadroniser::init()
    {
      pythia_->setLHAupPtr( static_cast<Pythia8::LHAup*>( lhaevt_.get() ) );
      if ( !pythia_->init() )
        FatalError( "Failed to initialise the Pythia8 core!\n\t"
                    "See the message above for more details." );
      return true;
    }

    void
    Pythia8Hadroniser::readString( const char* param )
    {
      if ( !pythia_->readString( param ) )
        FatalError( Form( "The Pythia8 core failed to parse the following setting:\n\t%s", param ) );
    }

    void
    Pythia8Hadroniser::setSeed( long long seed )
    {
#ifdef PYTHIA8
      if ( seed == -1ll ) {
        pythia_->settings.flag( "Random:setSeed", false );
        return;
      }
      pythia_->settings.flag( "Random:setSeed", true );
      pythia_->settings.parm( "Random:seed", seed );
#endif
    }

    bool
    Pythia8Hadroniser::decay( Event& ev, double& weight )
    {
      weight = 1.;
#ifndef PYTHIA8
      FatalError( "Pythia8 is not linked to this instance!" );
#else
      lhaevt_->feedEvent( ev );
      pythia_->next();
      lhaevt_->listEvent();
      pythia_->event.list();

      //updateEvent( ev, weight );
#endif
      return true;
    }

    bool
    Pythia8Hadroniser::hadronise( Event& ev, double& weight )
    {
      weight = 1.;
#ifndef PYTHIA8
      FatalError( "Pythia8 is not linked to this instance!" );
#else

      //===========================================================================================
      // launch the hadronisation / resonances decays
      //===========================================================================================

      if ( !launchPythia( ev ) )
        return false;

      updateEvent( ev, weight );

      ev.dump();
#endif
      return true;
    }

    bool
    Pythia8Hadroniser::launchPythia( Event& ev )
    {
      ev.num_hadronisation_trials = 0;
      while ( !pythia_->next() ) {
        //pythia_->event.list(true,true);
        //if ( pythia_->event.size() != num_py_parts ) break; //FIXME discards any pythia error!
        /*if ( proton_fragment ) {
          //std::cout << success << "attempt " << ev.num_hadronisation_trials << " / " << max_attempts_ << std::endl;
          pythia_->event.list(true,true);
        }*/
//        pythia_->stat();
        if ( ++ev.num_hadronisation_trials > max_attempts_ )
          return false;
          //exit(0);
      }
      return true;
    }

    void
    Pythia8Hadroniser::fragmentState( unsigned short idx, double xbj )
    {
      const Pythia8::Particle& remn = pythia_->event[idx];
      // specify the quark/diquark flavours
      //FIXME naive approach (weighted by e_q^2/1-e_dq^2); to be improved
      const double rnd = 1./RAND_MAX * rand();
      unsigned short pdg_q = 0, pdg_dq = 0;
      if      ( rnd < 1./9. ) { pdg_q = 1; pdg_dq = 2203; }
      else if ( rnd < 5./9. ) { pdg_q = 2; pdg_dq = 2101; }
      else                    { pdg_q = 2; pdg_dq = 2103; }
      // then assign the quark/diquark a 4-momentum
      const double px_x = remn.px(), px_y = remn.py(), px_z = remn.pz(), ex = remn.e();
      const double xdq = 1.-xbj;
      //const double m_q = pythia_->particleData.m0( pdg_q );
      //const double m_dq = pythia_->particleData.m0( pdg_dq );
      // fractional momenta of the two partons:
      // -> x * p_X for the quark
      // -> ( 1-x ) * p_X for the diquark
      Pythia8::Particle diquark( pdg_dq, 63, idx, 0, 0, 0, 0, 100+idx, px_x*xdq, px_y*xdq, px_z*xdq, ex*xdq );
      Pythia8::Particle   quark( pdg_q,  63, idx, 0, 0, 0, 100+idx, 0, px_x*xbj, px_y*xbj, px_z*xbj, ex*xbj );
      //diquark.e( diquark.eCalc() );
      //quark.e( quark.eCalc() );
      diquark.m( diquark.mCalc() );
      quark.m( quark.mCalc() );
      //std::cout << "> " << xdq << "|" << xbj << std::endl;
      const unsigned short id_dq = pythia_->event.append( diquark );
      const unsigned short id_q = pythia_->event.append( quark );
      // keep up with the particles parentage
      pythia_->event[idx].daughter1( id_dq );
      pythia_->event[idx].daughter2( id_q );
      //std::cout << "-->" << diquark.m() << "|" << quark.m() << std::endl;
      // set the quark/diquark to be hadronised through a string
      pythia_->event[idx].status( -15 );
    }

    void
    Pythia8Hadroniser::updateEvent( Event& ev, double& weight )
    {
      for ( unsigned short i = 1; i < pythia_->event.size(); ++i ) { // skip the central system
        const Pythia8::Particle& p = pythia_->event[i];
        if ( py_cg_corresp_.count( i ) > 0 ) {
          // the particle is already in the event content
          Particle& cg_part = ev.getById( py_cg_corresp_.at( i ) );
          if ( p.daughterList().size() == 0 )
            continue;
          if ( cg_part.role() != Particle::CentralSystem ) {
            cg_part.setStatus( Particle::Fragmented );
            continue;
          }
          weight *= p.particleDataEntry().pickChannel().bRatio();
          cg_part.setStatus( Particle::Resonance );
          continue;
        }
        // the particle was not yet included in the CepGen event
        const std::vector<int> mothers = p.motherList();
        if ( mothers.size() == 0 ) // isolated particle
          continue;

        Particle::Role role = Particle::CentralSystem;
        if ( py_cg_corresp_.count( mothers[0] ) > 0 ) {
          const Particle& moth = ev.getById( py_cg_corresp_.at( mothers[0] ) );
          if ( moth.role() == Particle::OutgoingBeam1
            || moth.role() == Particle::OutgoingBeam2 )
            role = moth.role();
        }

        Particle& op = ev.addParticle( role );
        cg_py_corresp_[op.id()] = i;
        py_cg_corresp_[i] = op.id();

        op.setPdgId( static_cast<ParticleCode>( abs( p.id() ) ), p.charge() );
        op.setStatus( p.isFinal()
          ? Particle::FinalState
          : Particle::Propagator
        );
        op.setMomentum( Particle::Momentum( p.px(), p.py(), p.pz(), p.e() ) );
        for ( const auto& moth : mothers ) {
          if ( moth != 0 && py_cg_corresp_.count( moth ) == 0 )
            FatalError( Form( "Particle with id=%d was not found in the event content!", moth ) );
          op.addMother( ev.getById( py_cg_corresp_.at( moth ) ) );
        }
      }
    }
  }

  LHAEvent::LHAEvent() :
    LHAup( 3 )
  {
    std::cout << __PRETTY_FUNCTION__ << std::endl;
    addProcess( 0, 1., 0., 1. );
  }

  void
  LHAEvent::feedEvent( const Event& ev )
  {
    setProcess( 0, 1., 100., Constants::alphaEM, Constants::alphaQCD );
    unsigned short parton1_id = 0, parton2_id = 0;
    for ( const auto& p : ev.particles() ) {
      int status = 0;
      unsigned short moth1 = 0, moth2 = 0;
      if ( p.mothers().size() > 0 ) {
        moth1 = *p.mothers().begin()+1;
        if ( p.mothers().size() > 1 ) {
          unsigned short tmp = *p.mothers().rbegin()+1;
          if ( tmp > moth1+1 ) {
            moth2 = moth1;
            moth1 = tmp;
          }
          else
            moth2 = tmp;
        }
      }
      const double mass2 = p.momentum().energy2()-p.momentum().p2(), mass = ( mass2 < 0 ) ? -sqrt( -mass2 ) : sqrt( mass2 );
      switch ( p.status() ) {
        case Particle::Unfragmented:
        case Particle::Undecayed:
          status = 1;
          break;
        case Particle::PrimordialIncoming:
          status = -9;
          break;
        case Particle::Incoming:
          status = 2;
          break;
        case Particle::FinalState:
        default:
          status = 3;
          break;
      }
      switch ( p.role() ) {
        case Particle::Parton1:
          parton1_id = sizePart();
          moth1 = moth2 = 0;
          addParticle( p.integerPdgId(), status, 0, 0, 0, 0, p.momentum().px(), p.momentum().py(), p.momentum().pz(), p.momentum().energy(), mass, 0., 0., 0. );
          break;
        case Particle::Parton2:
          parton2_id = sizePart();
          moth1 = moth2 = 0;
          addParticle( p.integerPdgId(), status, 0, 0, 0, 0, p.momentum().px(), p.momentum().py(), p.momentum().pz(), p.momentum().energy(), mass, 0., 0., 0. );
          break;
        case Particle::CentralSystem:
          addParticle( p.integerPdgId(), status, parton1_id, parton2_id, 0, 0, p.momentum().px(), p.momentum().py(), p.momentum().pz(), p.momentum().energy(), mass, 0., 0., 0. );
          continue;
        case Particle::IncomingBeam1:
          setBeamA( p.integerPdgId(), p.momentum().energy() );
          break;
        case Particle::IncomingBeam2:
          setBeamB( p.integerPdgId(), p.momentum().energy() );
          break;
        default:
          //break;
          continue;
      }
    }
  }

  bool
  LHAEvent::setInit()
  {
    return true;
  }

  bool
  LHAEvent::setEvent( int )
  {
    return true;
  }
}

#endif
