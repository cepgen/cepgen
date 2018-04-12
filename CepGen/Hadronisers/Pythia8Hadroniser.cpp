#include "CepGen/Hadronisers/Pythia8Hadroniser.h"

#include "CepGen/Parameters.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/PDG.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"

#include "CepGen/Event/Event.h"
#include "CepGen/Event/Particle.h"

#include "CepGen/Version.h"

namespace CepGen
{
  namespace Hadroniser
  {
    Pythia8::Vec4
    momToVec4( const Particle::Momentum& mom )
    {
//std::cout << mom << "|" << mom.mass() << "|" << Pythia8::Vec4( mom.px(), mom.py(), mom.pz(), mom.energy() ).mCalc() << std::endl;
      return Pythia8::Vec4( mom.px(), mom.py(), mom.pz(), mom.energy() );
    }

    Pythia8Hadroniser::Pythia8Hadroniser( const Parameters& params ) :
      GenericHadroniser( "pythia8" ), max_attempts_( params.hadroniser_max_trials ),
#ifdef PYTHIA8
      pythia_( new Pythia8::Pythia ), lhaevt_( new LHAEvent( &params ) ),
#endif
      full_evt_( false ), params_( &params )
    {
#ifdef PYTHIA8
      pythia_->setLHAupPtr( (Pythia8::LHAup*)lhaevt_.get() );
      pythia_->settings.parm( "Beams:idA", (short)params.kinematics.inpdg.first );
      pythia_->settings.parm( "Beams:idB", (short)params.kinematics.inpdg.second );
      // specify we will be using a LHA input
      pythia_->settings.mode( "Beams:frameType", 5 );
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
    Pythia8Hadroniser::init()
    {
      if ( pythia_->settings.flag( "ProcessLevel:all" ) != full_evt_ )
        pythia_->settings.flag( "ProcessLevel:all", full_evt_ );

      switch ( params_->kinematics.mode ) {
        case Kinematics::Mode::ElasticElastic:
          pythia_->settings.mode( "BeamRemnants:unresolvedHadron", 3 );
          break;
        case Kinematics::Mode::ElasticInelastic:
          pythia_->settings.mode( "BeamRemnants:unresolvedHadron", 1 );
          break;
        case Kinematics::Mode::InelasticElastic:
          pythia_->settings.mode( "BeamRemnants:unresolvedHadron", 2 );
          break;
        case Kinematics::Mode::InelasticInelastic: default:
          pythia_->settings.mode( "BeamRemnants:unresolvedHadron", 0 );
          break;
      }
//      pythia_->settings.mode( "BeamRemnants:remnantMode", 1 );

      if ( !pythia_->init() )
        throw CG_FATAL( "Pythia8Hadroniser" )
          << "Failed to initialise the Pythia8 core!\n\t"
          << "See the message above for more details.";

      return true;
    }

    void
    Pythia8Hadroniser::readString( const char* param )
    {
      if ( !pythia_->readString( param ) )
        throw CG_FATAL( "Pythia8Hadroniser" ) << "The Pythia8 core failed to parse the following setting:\n\t" << param;
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

full = true;//FIXME!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifndef PYTHIA8
      throw CG_FATAL( "Pythia8Hadroniser" ) << "Pythia8 is not linked to this instance!";
#else
      if ( !full && !pythia_->settings.flag( "ProcessLevel:resonanceDecays" ) )
        return true;

      if ( full != full_evt_ ) {
        full_evt_ = full;
        init();
      }

      //===========================================================================================
      // convert our event into a custom LHA format
      //===========================================================================================

      lhaevt_->feedEvent( ev, full );
lhaevt_->listEvent();
      //===========================================================================================
      // launch the hadronisation / resonances decays, and update the event accordingly
      //===========================================================================================

//lhaevt_->listEvent();
      unsigned short num_try = 0;
      while ( !pythia_->next() && num_try < ev.num_hadronisation_trials )
        num_try++;
//pythia_->process.list();
//pythia_->event.list();

      updateEvent( ev, weight, full );

#endif
      ev.dump();
      return true;
    }

    Particle&
    Pythia8Hadroniser::addParticle( Event& ev, const Pythia8::Particle& py_part, const Pythia8::Vec4& mom, unsigned short role, unsigned short offset ) const
    {
      Particle& op = ev.addParticle( (Particle::Role)role );
      op.setPdgId( static_cast<PDG>( abs( py_part.id() ) ), py_part.charge() );
      op.setStatus( py_part.isFinal()
        ? Particle::FinalState
        : Particle::Propagator
      );
      op.setMomentum( Particle::Momentum( mom.px(), mom.py(), mom.pz(), mom.e() ) );
      op.setMass( mom.mCalc() );
      lhaevt_->addCorresp( py_part.index()-offset, op.id() );
      return op;
    }

    void
    Pythia8Hadroniser::updateEvent( Event& ev, double& weight, bool full ) const
    {
      unsigned short offset = 0;

//pythia_->event.list();
      for ( unsigned short i = 1; i < pythia_->event.size(); ++i ) {
        if ( pythia_->event[i].status() != -12 )
          continue;
        // skip the incoming particles
        offset++;
      }

      for ( unsigned short i = 1+offset; i < pythia_->event.size(); ++i ) {
        const Pythia8::Particle& p = pythia_->event[i];
        const unsigned short cg_id = lhaevt_->cgPart( i-offset );
        if ( cg_id != LHAEvent::invalid_id ) { // particle already in the event
          Particle& cg_part = ev.getById( cg_id );
          if ( abs( p.id() ) != (unsigned short)cg_part.pdgId() )
            throw CG_FATAL( "Pythia8Hadroniser" ) << "Event list corruption detected for particle " << i << "!";
          if ( p.particleDataEntry().sizeChannels() == 0 )
            continue;
          weight *= p.particleDataEntry().pickChannel().bRatio();
          cg_part.setStatus( Particle::Resonance );
          continue;
        }
        // new particle to be added
        Particle::Role role = (Particle::Role)findRole( ev, p, offset );
        Pythia8::Vec4 mom = p.p();
        switch ( role ) {
          case Particle::OutgoingBeam1: // no break!
            ev.getByRole( Particle::OutgoingBeam1 )[0].setStatus( Particle::Fragmented );
            if ( abs( p.status() ) != 61 )
              break;
          case Particle::OutgoingBeam2: // no break!
            ev.getByRole( Particle::OutgoingBeam2 )[0].setStatus( Particle::Fragmented );
            if ( abs( p.status() ) != 61 )
              break;
          default:
            break;
        }
        // found the role ; now we can add the particle
        Particle& cg_part = addParticle( ev, p, mom, (unsigned short)role, offset );
        for ( const auto& moth_id : p.motherList() ) {
          if ( moth_id <= offset )
            continue;
          const unsigned short moth_cg_id = lhaevt_->cgPart( moth_id-offset );
          if ( moth_cg_id != LHAEvent::invalid_id ) {
            Particle& moth = ev.getById( moth_cg_id );
            cg_part.addMother( moth );
          }
          else {
            Particle& cg_moth = addParticle( ev, pythia_->event[moth_id], mom, (unsigned short)role, offset );
            cg_part.addMother( cg_moth );
          }
        }
      }
    }

    unsigned short
    Pythia8Hadroniser::findRole( const Event& ev, const Pythia8::Particle& p, unsigned short offset ) const
    {
      for ( const auto& par_id : p.motherList() ) {
        if ( offset > 0 && par_id == 1 )
          return (unsigned short)Particle::OutgoingBeam1;
        if ( offset > 0 && par_id == 2 )
          return (unsigned short)Particle::OutgoingBeam2;
        const unsigned short par_cg_id = lhaevt_->cgPart( par_id-offset );
        if ( par_cg_id != LHAEvent::invalid_id )
          return (unsigned short)ev.getConstById( par_cg_id ).role();
        return findRole( ev, pythia_->event[par_id], offset );
      }
      return (unsigned short)Particle::UnknownRole;
    }
  }

  //================================================================================================
  // Custom LHA event definition
  //================================================================================================

#ifdef PYTHIA8
  const double LHAEvent::mp_ = ParticleProperties::mass( PDG::Proton );
  const double LHAEvent::mp2_ = LHAEvent::mp_*LHAEvent::mp_;

  LHAEvent::LHAEvent( const Parameters* params ) :
    LHAup( 3 ), params_( params )
  {
    addProcess( 0, 10., 0., 100. );
    if ( params_ ) {
      setBeamA( (short)params_->kinematics.inpdg.first, params_->kinematics.inp.first );
      setBeamB( (short)params_->kinematics.inpdg.second, params_->kinematics.inp.second );
    }
  }

  void
  LHAEvent::setCrossSection( int id, double xsec, double xsec_err )
  {
    setXSec( id, xsec );
    setXErr( id, xsec_err );
    //listInit();
  }

  void
  LHAEvent::feedEvent( const Event& ev, bool full )
  {
    const double scale = ev.getOneByRole( Particle::Intermediate ).mass();
    setProcess( 0, 1., scale, Constants::alphaEM, Constants::alphaQCD );

    Pythia8::Vec4 mom_iq1, mom_iq2;
    unsigned short quark1_id = 0, quark2_id = 0, quark1_pdgid = 0, quark2_pdgid = 0;
    double x1 = 0., x2 = 0.;

    const Particle& part1 = ev.getOneByRole( Particle::Parton1 ), &part2 = ev.getOneByRole( Particle::Parton2 );

    if ( full ) {
      const Particle& op1 = ev.getOneByRole( Particle::OutgoingBeam1 ), &op2 = ev.getOneByRole( Particle::OutgoingBeam2 );
      const Particle& ip1 = ev.getOneByRole( Particle::IncomingBeam1 ), &ip2 = ev.getOneByRole( Particle::IncomingBeam2 );
      const double q2_1 = -part1.momentum().mass2(), q2_2 = -part2.momentum().mass2();
      x1 = q2_1/( q2_1+op1.mass2()-mp2_ );
      x2 = q2_2/( q2_2+op2.mass2()-mp2_ );

      quark1_pdgid = quark2_pdgid = 2; //FIXME

      setIdX( op1.integerPdgId(), op2.integerPdgId(), x1, x2 ); //FIXME initiator = photon or proton?
      //setIdX( quark1_pdgid, quark2_pdgid, x1, x2 );

      //===========================================================================================
      // incoming valence quarks
      //===========================================================================================

      mom_iq1 = Pythia8::Vec4( 0., 0., x1*ip1.momentum().pz(), x1*ip1.energy() );
      quark1_id = sizePart();
      addCorresp( quark1_id, op1.id() );
      addParticle( quark1_pdgid, -1, 0, 0, 501, 0, mom_iq1.px(), mom_iq1.py(), mom_iq1.pz(), mom_iq1.e(), mom_iq1.mCalc(), 0., 1. );

      mom_iq2 = Pythia8::Vec4( 0., 0., x2*ip2.momentum().pz(), x2*ip2.energy() );
      quark2_id = sizePart();
      addCorresp( quark2_id, op2.id() );
      addParticle( quark2_pdgid, -1, 0, 0, 502, 0, mom_iq2.px(), mom_iq2.py(), mom_iq2.pz(), mom_iq2.e(), mom_iq2.mCalc(), 0., 1. );
      /*const Pythia8::Vec4 mom_dq1( Hadroniser::momToVec4( ip1.momentum() )-mom_iq1 );
      addParticle( 2101, -1, 0, 0, 0, 501, mom_dq1.px(), mom_dq1.py(), mom_dq1.pz(), mom_dq1.e(), mom_dq1.mCalc(), 0., 1. );
      const Pythia8::Vec4 mom_dq2( Hadroniser::momToVec4( ip2.momentum() )-mom_iq2 );
      addParticle( 2101, -1, 0, 0, 0, 502, mom_dq2.px(), mom_dq2.py(), mom_dq2.pz(), mom_dq2.e(), mom_dq2.mCalc(), 0., 1. );*/
    }
//    std::cout << x1 << "|" << x2 << std::endl;

    //=============================================================================================
    // incoming partons
    //=============================================================================================

    const Pythia8::Vec4 mom_part1( Hadroniser::momToVec4( part1.momentum() ) );
    addCorresp( sizePart(), part1.id() );
    addParticle( part1.integerPdgId(), 2, quark1_id, quark1_id+2, 0, 0, mom_part1.px(), mom_part1.py(), mom_part1.pz(), mom_part1.e(), mom_part1.mCalc(), 0., 0. );

    const Pythia8::Vec4 mom_part2( Hadroniser::momToVec4( part2.momentum() ) );
    addCorresp( sizePart(), part2.id() );
    addParticle( part2.integerPdgId(), 2, quark2_id, quark2_id+2, 0, 0, mom_part2.px(), mom_part2.py(), mom_part2.pz(), mom_part2.e(), mom_part2.mCalc(), 0., 0. );

/*    if ( full ) {
      //===========================================================================================
      // outgoing valence quarks
      //===========================================================================================

      const Pythia8::Vec4 mom_oq1 = mom_iq1-mom_part1;
      //addParticle( quark1_pdgid, 1, quark1_id, pyPart( part1.id() ), 501, 0, mom_oq1.px(), mom_oq1.py(), mom_oq1.pz(), mom_oq1.e(), mom_oq1.mCalc() );
      addParticle( quark1_pdgid, 1, quark1_id, 0, 501, 0, mom_oq1.px(), mom_oq1.py(), mom_oq1.pz(), mom_oq1.e(), mom_oq1.mCalc(), 0., 1. );

      const Pythia8::Vec4 mom_oq2 = mom_iq2-mom_part2;
      //addParticle( quark2_pdgid, 1, quark2_id, pyPart( part2.id() ), 502, 0, mom_oq2.px(), mom_oq2.py(), mom_oq2.pz(), mom_oq2.e(), mom_oq2.mCalc() );
      addParticle( quark2_pdgid, 1, quark2_id, 0, 502, 0, mom_oq2.px(), mom_oq2.py(), mom_oq2.pz(), mom_oq2.e(), mom_oq2.mCalc(), 0., 1. );
//std::cout << "MX=" << op1.momentum().mass() << "|MY=" << op2.momentum().mass() << std::endl;
    }*/

    //=============================================================================================
    // central system
    //=============================================================================================

    for ( const auto& p : ev.getByRole( Particle::CentralSystem ) ) {
      addCorresp( sizePart(), p.id() );
      Pythia8::Vec4 mom_part( p.momentum().px(), p.momentum().py(), p.momentum().pz(), p.momentum().energy() );
      const auto mothers = p.mothers();
      //unsigned short moth1_id = 1, moth2_id = 2;
      unsigned short moth1_id = 0, moth2_id = 0;
      if ( mothers.size() > 0 ) {
        const unsigned short moth1_cg_id = *mothers.begin();
        moth1_id = pyPart( moth1_cg_id );
        if ( moth1_id == invalid_id ) {
          const Particle& moth = ev.getConstById( moth1_cg_id );
          if ( moth.mothers().size() > 0 )
            moth1_id = pyPart( *moth.mothers().begin() );
          if ( moth.mothers().size() > 1 )
            moth2_id = pyPart( *moth.mothers().rbegin() );
        }
        if ( mothers.size() > 1 ) {
          const unsigned short moth2_cg_id = *mothers.rbegin();
          moth2_id = pyPart( moth2_cg_id );
          if ( moth2_id == invalid_id ) {
            const Particle& moth = ev.getConstById( moth2_cg_id );
            moth.dump();
            moth2_id = pyPart( *moth.mothers().rbegin() );
          }
        }
      }
      addParticle( p.integerPdgId(), 1, moth1_id, moth2_id, 0, 0, mom_part.px(), mom_part.py(), mom_part.pz(), mom_part.e(), mom_part.mCalc(), 0., 0., 0. );
    }
    if ( full )
      //setPdf( quark1_pdgid, quark2_pdgid, x1, x2, scale, x1*( 1-x1 ), x2*( 1-x2 ), true );
      setPdf( quark1_pdgid, quark2_pdgid, x1, x2, scale, 0., 0., false );
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
    CG_INFO( "LHAEvent" ) << oss.str();
  }
#endif
}

