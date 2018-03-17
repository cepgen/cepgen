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
      GenericHadroniser( "pythia8" ), max_attempts_( params.hadroniser_max_trials ),
#ifdef PYTHIA8
      pythia_( new Pythia8::Pythia ), lhaevt_( new LHAEvent( &params ) ),
#endif
      full_evt_( false ), params_( &params )
    {
#ifdef PYTHIA8
      pythia_->setLHAupPtr( (Pythia8::LHAup*)lhaevt_.get() );
      pythia_->settings.parm( "Beams:idA", params.kinematics.inpdg.first );
      pythia_->settings.parm( "Beams:idB", params.kinematics.inpdg.second );
      pythia_->settings.parm( "Beams:eCM", params.kinematics.sqrtS() );
#endif
      for ( const auto& pdgid : params.kinematics.minimum_final_state )
        min_ids_.emplace_back( (unsigned short)pdgid );
    }

    Pythia8Hadroniser::~Pythia8Hadroniser()
    {
      pythia_->stat();
    }

    bool
    Pythia8Hadroniser::init( bool enable_all_processes )
    {
      //enable_all_processes = true;//FIXME FIXME
      if ( pythia_->settings.flag( "ProcessLevel:all" ) != enable_all_processes )
        pythia_->settings.flag( "ProcessLevel:all", enable_all_processes );

      if ( !pythia_->init() )
        FatalError( "Failed to initialise the Pythia8 core!\n\t"
                    "See the message above for more details." );

      full_evt_ = enable_all_processes;
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
      if ( full && !full_evt_ )
        init( true );

      //===========================================================================================
      // convert our event into a custom LHA format
      //===========================================================================================
//      std::cout << "bbb" << std::endl;
      lhaevt_->feedEvent( ev, full );
//      lhaevt_->listEvent();

      //===========================================================================================
      // launch the hadronisation / resonances decays, and update the event accordingly
      //===========================================================================================

      unsigned short num_try = 0;
      while ( !pythia_->next() && num_try < ev.num_hadronisation_trials )
        num_try++;

//      pythia_->event.list();
      updateEvent( pythia_.get(), lhaevt_.get(), ev, weight );
      //lhaevt_->listEvent();
//      if ( pythia_->info.nMPI() > 0 )
//        std::cout << "#MPI:" << pythia_->info.nMPI() << std::endl;

#endif
      return true;
    }

    Particle&
    Pythia8Hadroniser::addParticle( LHAEvent* lhaevt, Event& ev, const Pythia8::Particle& py_part, unsigned short role, unsigned short offset ) const
    {
      Particle& op = ev.addParticle( (Particle::Role)role );
      op.setPdgId( static_cast<ParticleCode>( abs( py_part.id() ) ), py_part.charge() );
      op.setStatus( py_part.isFinal()
        ? Particle::FinalState
        : Particle::Propagator
      );
      op.setMomentum( Particle::Momentum( py_part.px(), py_part.py(), py_part.pz(), py_part.e() ) );
      op.setMass( py_part.m() );
      lhaevt->addCorresp( py_part.index()-offset, op.id() );
      return op;
    }

    void
    Pythia8Hadroniser::updateEvent( const Pythia8::Pythia* pythia, LHAEvent* lhaevt, Event& ev, double& weight ) const
    {
      unsigned short offset = 0;
      for ( unsigned short i = 1; i < pythia->event.size(); ++i ) {
        const Pythia8::Particle& p = pythia->event[i];
        if ( p.status() == -12 ) { // skip the incoming particles
          offset++;
          continue;
        }
        const unsigned short cg_id = lhaevt->cgPart( i-offset );
        if ( cg_id != LHAEvent::invalid_id ) { // particle already in the event
          Particle& cg_part = ev.getById( cg_id );
          if ( abs( p.id() ) != (unsigned short)cg_part.pdgId() )
            FatalError( Form( "Event list corruption detected for particle %d", i ) );
          if ( p.particleDataEntry().sizeChannels() == 0 )
            continue;
          weight *= p.particleDataEntry().pickChannel().bRatio();
          cg_part.setStatus( Particle::Resonance );
          continue;
        }
        // new particle to be added
        Particle::Role role = (Particle::Role)findRole( pythia, lhaevt, ev, p, offset );
        if ( role == Particle::OutgoingBeam1 )
          ev.getByRole( Particle::OutgoingBeam1 )[0].setStatus( Particle::Fragmented );
        if ( role == Particle::OutgoingBeam2 )
          ev.getByRole( Particle::OutgoingBeam2 )[0].setStatus( Particle::Fragmented );
        // found the role ; now we can add the particle
        Particle& cg_part = addParticle( lhaevt, ev, p, (unsigned short)role, offset );
        for ( const auto& moth_id : p.motherList() ) {
          if ( moth_id <= offset )
            continue;
          const unsigned short moth_cg_id = lhaevt->cgPart( moth_id-offset );
          if ( moth_cg_id != LHAEvent::invalid_id ) {
            Particle& moth = ev.getById( moth_cg_id );
            cg_part.addMother( moth );
          }
          else {
            Particle& cg_moth = addParticle( lhaevt, ev, pythia->event[moth_id], (unsigned short)role, offset );
            cg_part.addMother( cg_moth );
          }
        }
      }
//      lhaevt_->dumpCorresp();
//      ev.dump();
//      exit(0);
    }

    unsigned short
    Pythia8Hadroniser::findRole( const Pythia8::Pythia* pythia, const LHAEvent* lhaevt, const Event& ev, const Pythia8::Particle& p, unsigned short offset ) const
    {
      for ( const auto& par_id : p.motherList() ) {
        if ( offset > 0 && par_id == 1 )
          return (unsigned short)Particle::OutgoingBeam1;
        if ( offset > 0 && par_id == 2 )
          return (unsigned short)Particle::OutgoingBeam2;
        const unsigned short par_cg_id = lhaevt->cgPart( par_id-offset );
        if ( par_cg_id != LHAEvent::invalid_id )
          return (unsigned short)ev.getConstById( par_cg_id ).role();
        return findRole( pythia, lhaevt, ev, pythia->event[par_id], offset );
      }
      return (unsigned short)Particle::UnknownRole;
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
//    setIdX( 2212, 2212, 0., 0. );
//    setIdX( 22, 22, 0.5, 0.5 );
    if ( params_ ) {
      setBeamA( (short)params_->kinematics.inpdg.first, params_->kinematics.inp.first );
      setBeamB( (short)params_->kinematics.inpdg.second, params_->kinematics.inp.second );
    }
  }

  void
  LHAEvent::setCrossSection( int id, double xsec, double xsec_err )
  {
//    addProcess( id, xsec, xsec_err );
    setXSec( id, xsec );
    setXErr( id, xsec_err );
    listInit();
  }

  void
  LHAEvent::feedEvent( const Event& ev, bool full )
  {
    const double scale = ev.getOneByRole( Particle::Intermediate ).mass();
    setProcess( 0, 1., scale, Constants::alphaEM, Constants::alphaQCD );

    double x1 = 0., x2 = 0.;
    if ( full ) {
      const Particle& ip1 = ev.getOneByRole( Particle::IncomingBeam1 ), &ip2 = ev.getOneByRole( Particle::IncomingBeam2 );
      const Particle& op1 = ev.getOneByRole( Particle::OutgoingBeam1 ), &op2 = ev.getOneByRole( Particle::OutgoingBeam2 );
      const double q2_1 = -( ip1.momentum()-op1.momentum() ).mass2(), x1 = q2_1/( q2_1+op1.mass2()-mp2_ );
      const double q2_2 = -( ip2.momentum()-op2.momentum() ).mass2(), x2 = q2_2/( q2_2+op2.mass2()-mp2_ );
      //const double x1 = ( ip1.momentum()-op1.momentum() ).energy()/ip1.energy(), x2 = ( ip2.momentum()-op2.momentum() ).energy()/ip2.energy();
      setIdX( ip1.integerPdgId(), ip2.integerPdgId(), x1, x2 );
    }

    unsigned short parton1_id = 0, parton2_id = 0, parton1_pdgid = 0, parton2_pdgid = 0;
    for ( const auto& p : ev.particles() ) {
      switch ( p.role() ) {
        case Particle::Parton1:
        case Particle::Parton2: {
          if ( p.role() == Particle::Parton1 ) {
            parton1_id = sizePart();
            parton1_pdgid = p.integerPdgId();
          }
          else {
            parton2_id = sizePart();
            parton2_pdgid = p.integerPdgId();
          }
          py_cg_corresp_.emplace_back( sizePart(), p.id() );
          if ( full )
            addParticle( p.integerPdgId(), -2, 0, 0, 0, 0, p.momentum().px(), p.momentum().py(), p.momentum().pz(), p.momentum().energy(), p.mass(), 0., 0., -p.momentum().mass() );
          else
            addParticle( p.integerPdgId(), -2, 0, 0, 0, 0, p.momentum().px(), p.momentum().py(), p.momentum().pz(), p.momentum().energy(), p.momentum().mass(), 0., 0., 0. );
        } break;
        case Particle::CentralSystem: {
          py_cg_corresp_.emplace_back( sizePart(), p.id() );
          addParticle( p.integerPdgId(), 1, parton1_id, parton2_id, 0, 0, p.momentum().px(), p.momentum().py(), p.momentum().pz(), p.momentum().energy(), p.momentum().mass(), 0., 0., 0. );
          continue;
        } break;
        default:
          continue;
      }
    }
    if ( full )
      setPdf( parton1_pdgid, parton2_pdgid, x1, x2, scale, 0., 0., true );
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
}

#endif
