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
      lhaevt_( new LHAEvent( &params ) )
    {
#ifdef PYTHIA8
      pythia_.reset( new Pythia8::Pythia );
      pythia_->settings.parm( "Beams:idA", params.kinematics.inpdg.first );
      pythia_->settings.parm( "Beams:idB", params.kinematics.inpdg.second );
      pythia_->settings.parm( "Beams:eCM", params.kinematics.sqrtS() );
      pythia_->setLHAupPtr( static_cast<Pythia8::LHAup*>( lhaevt_.get() ) );
#endif
      for ( const auto& pdgid : params.kinematics.minimum_final_state )
        min_ids_.emplace_back( (unsigned short)pdgid );
    }

    Pythia8Hadroniser::~Pythia8Hadroniser()
    {}

    bool
    Pythia8Hadroniser::init()
    {
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
      pythia_->settings.mode( "Random:seed", seed );
#endif
    }

    void
    Pythia8Hadroniser::setCrossSection( double xsec, double xsec_err )
    {
      lhaevt_->setCrossSection( 0, xsec, xsec_err );
    }

    bool
    Pythia8Hadroniser::run( Event& ev, double& weight, bool full )
    {
      weight = 1.;
#ifndef PYTHIA8
      FatalError( "Pythia8 is not linked to this instance!" );
#else
      //===========================================================================================
      // convert our event into a custom LHA format
      //===========================================================================================

      lhaevt_->feedEvent( ev, full );
//      lhaevt_->feedEvent( ev, false );

      //===========================================================================================
      // launch the hadronisation / resonances decays
      //===========================================================================================

      pythia_->next();
      //lhaevt_->listEvent();
      pythia_->event.list();
      if ( pythia_->info.nMPI() > 0 )
        std::cout << "#MPI:" << pythia_->info.nMPI() << std::endl;

      //===========================================================================================
      // update our event
      //===========================================================================================

      updateEvent( ev, weight );
#endif
      return true;
    }

    bool
    Pythia8Hadroniser::launchPythia( Event& ev )
    {
      ev.num_hadronisation_trials = 0;
      while ( !pythia_->next() ) {
        if ( ++ev.num_hadronisation_trials > max_attempts_ )
          return false;
        //pythia_->event.list(true,true);
        //if ( pythia_->event.size() == num_py_parts )
        //  continue; //FIXME discards any pythia error!
        /*if ( proton_fragment ) {
          //std::cout << success << "attempt " << ev.num_hadronisation_trials << " / " << max_attempts_ << std::endl;
          pythia_->event.list(true,true);
        }*/
//        pythia_->stat();
          //exit(0);
      }
      return true;
    }

    void
    Pythia8Hadroniser::updateEvent( Event& ev, double& weight )
    {
      ev.dump();
      lhaevt_->dumpCorresp();
      for ( unsigned short i = 1; i < pythia_->event.size(); ++i ) { // skip the central system
        const Pythia8::Particle& p = pythia_->event[i];
        const unsigned short cg_id = lhaevt_->cgPart( i );
        if ( cg_id != LHAEvent::invalid_id ) {
          // the particle is already in the event content
          Particle& cg_part = ev.getById( cg_id );
          cg_part.dump();
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
        const unsigned short moth_id = lhaevt_->cgPart( mothers[0] );
        if ( moth_id != LHAEvent::invalid_id ) {
          const Particle& moth = ev.getById( moth_id );
          if ( moth.role() == Particle::OutgoingBeam1
            || moth.role() == Particle::OutgoingBeam2 )
            role = moth.role();
        }

        Particle& op = ev.addParticle( role );
        lhaevt_->addCorresp( i, op.id() );

        op.setPdgId( static_cast<ParticleCode>( abs( p.id() ) ), p.charge() );
        op.setStatus( p.isFinal()
          ? Particle::FinalState
          : Particle::Propagator
        );
        op.setMomentum( Particle::Momentum( p.px(), p.py(), p.pz(), p.e() ) );
        for ( const auto& moth : mothers ) {
          const unsigned short moth_id = lhaevt_->cgPart( moth );
          if ( moth != 0 && moth_id == LHAEvent::invalid_id )
            FatalError( Form( "Particle with id=%d was not found in the event content!", moth ) );
          op.addMother( ev.getById( moth_id ) );
        }
      }
      ev.dump();
    }
  }

  //================================================================================================
  // Custom LHA event definition
  //================================================================================================

  const double LHAEvent::mp_ = ParticleProperties::mass( Proton );
  const double LHAEvent::mp2_ = LHAEvent::mp_*LHAEvent::mp_;

  LHAEvent::LHAEvent( const Parameters* params ) :
    LHAup( 3 ), params_( params )
  {
    addProcess( 0, 10., 0., 100. );
    setBeamA( (short)params_->kinematics.inpdg.first, params_->kinematics.inp.first );
    setBeamB( (short)params_->kinematics.inpdg.second, params_->kinematics.inp.second );
  }

  void
  LHAEvent::setCrossSection( int id, double xsec, double xsec_err )
  {
    setXSec( id, xsec );
    setXErr( id, xsec_err );
    listInit();
  }

  void
  LHAEvent::feedEvent( const Event& ev, bool full )
  {
    setProcess( 0, 1., ev.getOneByRole( Particle::Intermediate ).mass(), Constants::alphaEM, Constants::alphaQCD );
    unsigned short parton1_id = 0, parton2_id = 0;
    double x1 = 0., x2 = 0.;
    int pdg1 = 2212, pdg2 = 2212;
    for ( const auto& p : ev.particles() ) {
      int status = 0;
      switch ( p.status() ) {
        case Particle::Unfragmented:
          status = ( full ) ? -2 : 1;
          break;
        case Particle::Undecayed:
          status = 1;
          break;
        case Particle::Incoming:
          status = -1;
          break;
        case Particle::FinalState:
        default:
          break;
      }
      switch ( p.role() ) {
        case Particle::Parton1: {
          parton1_id = sizePart();
          py_cg_corresp_.emplace_back( sizePart(), p.id() );
          /*unsigned short moth = ( p.mothers().size() == 0 ) ? 0 : pyPart( *p.mothers().begin() );
          if ( moth == invalid_id )
            moth = 0;
          addParticle( p.integerPdgId(), status, moth, 0, 0, 0, p.momentum().px(), p.momentum().py(), p.momentum().pz(), p.momentum().energy(), mass, 0., 0., 0. );*/
          addParticle( p.integerPdgId(), status, 0, 0, 0, 0, p.momentum().px(), p.momentum().py(), p.momentum().pz(), p.momentum().energy(), p.momentum().mass(), 0., 0., 0. );
        } break;
        case Particle::Parton2: {
          parton2_id = sizePart();
          py_cg_corresp_.emplace_back( sizePart(), p.id() );
          /*unsigned short moth = ( p.mothers().size() == 0 ) ? 0 : pyPart( *p.mothers().begin() );
          if ( moth == invalid_id )
            moth = 0;
          addParticle( p.integerPdgId(), status, moth, 0, 0, 0, p.momentum().px(), p.momentum().py(), p.momentum().pz(), p.momentum().energy(), mass, 0., 0., 0. );*/
          addParticle( p.integerPdgId(), status, 0, 0, 0, 0, p.momentum().px(), p.momentum().py(), p.momentum().pz(), p.momentum().energy(), p.momentum().mass(), 0., 0., 0. );
        } break;
        case Particle::CentralSystem:
          py_cg_corresp_.emplace_back( sizePart(), p.id() );
          addParticle( p.integerPdgId(), status, parton1_id, parton2_id, 0, 0, p.momentum().px(), p.momentum().py(), p.momentum().pz(), p.momentum().energy(), p.momentum().mass(), 0., 0., 0. );
          continue;
        case Particle::OutgoingBeam1: {
          if ( !full )
            continue;
          const Particle& parent = ev.getConstById( *p.mothers().begin() );
          const double q2 = -( parent.momentum()-p.momentum() ).mass2();
          pdg1 = parent.integerPdgId();
          x1 = q2/( q2+p.mass2()-mp2_ );
          //fragmentState( p, x1 );
        } break;
        case Particle::OutgoingBeam2: {
          if ( !full )
            continue;
          const Particle& parent = ev.getConstById( *p.mothers().begin() );
          const double q2 = -( parent.momentum()-p.momentum() ).mass2();
          pdg2 = parent.integerPdgId();
          x2 = q2/( q2+p.mass2()-mp2_ );
          //fragmentState( p, x2 );
        } break;
        /*case Particle::IncomingBeam1:
          setBeamA( p.integerPdgId(), p.momentum().energy() );
          py_cg_corresp_.emplace_back( sizePart(), p.id() );
          addParticle( p.integerPdgId(), -9, 0, 0, 0, 0, p.momentum().px(), p.momentum().py(), p.momentum().pz(), p.momentum().energy(), mass, 0., 0., 0. );
          break;
        case Particle::IncomingBeam2:
          setBeamB( p.integerPdgId(), p.momentum().energy() );
          py_cg_corresp_.emplace_back( sizePart(), p.id() );
          addParticle( p.integerPdgId(), -9, 0, 0, 0, 0, p.momentum().px(), p.momentum().py(), p.momentum().pz(), p.momentum().energy(), mass, 0., 0., 0. );
          break;*/
        default:
          //break;
          continue;
      }
    }
    setIdX( pdg1, pdg2, x1, x2 );
//    listEvent();
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

  void
  LHAEvent::setProcess( int id, double xsec, double q2_scale, double alpha_qed, double alpha_qcd )
  {
    LHAup::setProcess( id, xsec, q2_scale, alpha_qed, alpha_qcd );
    py_cg_corresp_.clear();
  }

  unsigned short
  LHAEvent::cgPart( unsigned short py_id ) const
  {
    for ( const auto& py_cg : py_cg_corresp_ )
      if ( py_cg.first == py_id )
        return py_cg.second;
    return invalid_id;
  }

  unsigned short
  LHAEvent::pyPart( unsigned short cg_id ) const
  {
    for ( const auto& py_cg : py_cg_corresp_ )
      if ( py_cg.second == cg_id )
        return py_cg.first;
    return invalid_id;
  }

  void
  LHAEvent::addCorresp( unsigned short py_id, unsigned short cg_id )
  {
    py_cg_corresp_.emplace_back( py_id, cg_id );
  }

  void
  LHAEvent::dumpCorresp() const
  {
    std::ostringstream oss;
    oss << "List of Pythia <-> CepGen particle ids correspondance";
    for ( const auto& py_cg : py_cg_corresp_ )
      oss << "\n\t" << py_cg.first << " <-> " << py_cg.second;
    Information( oss.str() );
  }

  void
  LHAEvent::fragmentState( const Particle& p_x, double xbj )
  {
    unsigned short pdg_q = 0, pdg_dq = 0;
    {
      // specify the quark/diquark flavours
      //FIXME naive approach (weighted by e_q^2/1-e_dq^2); to be improved
      const double rnd = 1./RAND_MAX * rand();
      if      ( rnd < 1./9. ) { pdg_q = 1; pdg_dq = 2203; }
      else if ( rnd < 5./9. ) { pdg_q = 2; pdg_dq = 2101; }
      else                    { pdg_q = 2; pdg_dq = 2103; }
    }
    // then assign the quark/diquark a 4-momentum
    const Particle::Momentum p_q = p_x.momentum()*xbj, p_dq = p_x.momentum()-p_q;
    // fractional momenta of the two partons:
    // -> x * p_X for the quark
    // -> ( 1-x ) * p_X for the diquark
//    std::cout << p_x.momentum() << "|" << p_q << "|" << p_dq << std::endl;
    py_cg_corresp_.emplace_back( sizePart(), p_x.id() );
    addParticle( pdg_dq, 1, 0, 0, 0, 100+p_x.id(), p_dq.px(), p_dq.py(), p_dq.pz(), p_dq.energy(), p_dq.mass() );
    py_cg_corresp_.emplace_back( sizePart(), p_x.id() );
    addParticle( pdg_q, 1, 0, 0, 100+p_x.id(), 0, p_q.px(), p_q.py(), p_q.pz(), p_q.energy(), p_q.mass() );
  }
}

#endif
