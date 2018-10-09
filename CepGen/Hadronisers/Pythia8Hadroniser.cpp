#include "CepGen/Hadronisers/Pythia8Hadroniser.h"

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
#ifdef PYTHIA8
    Pythia8::Vec4
    momToVec4( const Particle::Momentum& mom )
    {
      return Pythia8::Vec4( mom.px(), mom.py(), mom.pz(), mom.energy() );
    }
#endif

    Pythia8Hadroniser::Pythia8Hadroniser( const Parameters& params, const ParametersList& plist ) :
      GenericHadroniser( "pythia8", plist ),
#ifdef PYTHIA8
      pythia_( new Pythia8::Pythia ), lhaevt_( new LHAEvent( &params ) ),
#endif
      full_evt_( false ), offset_( 0 ), first_evt_( true ), params_( &params )
    {
#ifdef PYTHIA8
      pythia_->setLHAupPtr( (Pythia8::LHAup*)lhaevt_.get() );
      pythia_->settings.parm( "Beams:idA", (short)params.kinematics.incoming_beams.first.pdg );
      pythia_->settings.parm( "Beams:idB", (short)params.kinematics.incoming_beams.second.pdg );
      // specify we will be using a LHA input
      pythia_->settings.mode( "Beams:frameType", 5 );
      pythia_->settings.parm( "Beams:eCM", params.kinematics.sqrtS() );
#endif
      for ( const auto& pdgid : params.kinematics.minimum_final_state )
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
      lhaevt_->setCrossSection( 0, xsec, xsec_err );
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

      lhaevt_->feedEvent( ev, full, params_->kinematics.mode );
      //if ( full ) lhaevt_->listEvent();

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

      updateEvent( ev, weight, full );
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
      op.setPdgId( static_cast<PDG>( abs( py_part.id() ) ), py_part.charge() );
      op.setStatus( py_part.isFinal()
        ? Particle::Status::FinalState
        : Particle::Status::Propagator );
      op.setMomentum( Particle::Momentum( mom.px(), mom.py(), mom.pz(), mom.e() ) );
      op.setMass( mom.mCalc() );
      lhaevt_->addCorresp( py_part.index()-offset_, op.id() );
      return op;
    }

    void
    Pythia8Hadroniser::updateEvent( Event& ev, double& weight, bool full ) const
    {
      for ( unsigned short i = 1+offset_; i < pythia_->event.size(); ++i ) {
        const Pythia8::Particle& p = pythia_->event[i];
        const unsigned short cg_id = lhaevt_->cepgenId( i-offset_ );
        if ( cg_id != LHAEvent::invalid_id ) {
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
            lhaevt_->listEvent();
            CG_INFO( "Pythia8Hadroniser:update" ) << "Pythia event content:";
            pythia_->event.list();
            CG_INFO( "Pythia8Hadroniser:update" ) << "CepGen event content:";
            ev.dump();
            CG_INFO( "Pythia8Hadroniser:update" ) << "Correspondence:";
            lhaevt_->dumpCorresp();

            throw CG_FATAL( "Pythia8Hadroniser:update" )
              << "Event list corruption detected for (Pythia/CepGen) particle " << i << "/" << cg_id << ":\n\t"
              << "should be " << abs( p.id() ) << ", "
              << "got " << cg_part.integerPdgId() << "!";
          }
          //--- resonance decayed; apply branching ratio for this decay
          if ( p.particleDataEntry().sizeChannels() > 0 ) {
            weight *= p.particleDataEntry().pickChannel().bRatio();
            cg_part.setStatus( Particle::Status::Resonance );
          }
        }
        else {
          //----- new particle to be added
          const unsigned short role = findRole( ev, p );
          switch ( (Particle::Role)role ) {
            default: break;
            case Particle::OutgoingBeam1: {
              ev.getByRole( Particle::OutgoingBeam1 )[0].setStatus( Particle::Status::Fragmented );
              if ( abs( p.status() ) != 61 )
                break;
            } // no break!
            case Particle::OutgoingBeam2: {
              ev.getByRole( Particle::OutgoingBeam2 )[0].setStatus( Particle::Status::Fragmented );
              if ( abs( p.status() ) != 61 )
                break;
            } // no break!
          }
          // found the role ; now we can add the particle
          Particle& cg_part = addParticle( ev, p, p.p(), role );
          for ( const auto& moth_id : p.motherList() ) {
            if ( moth_id <= offset_ )
              continue;
            const unsigned short moth_cg_id = lhaevt_->cepgenId( moth_id-offset_ );
            if ( moth_cg_id != LHAEvent::invalid_id )
              cg_part.addMother( ev[moth_cg_id] );
            else
              cg_part.addMother( addParticle( ev, pythia_->event[moth_id], p.p(), role ) );
            if ( !p.isFinal() ) {
              if ( p.isResonance() || p.daughterList().size() > 0 )
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
        const unsigned short par_cg_id = lhaevt_->cepgenId( par_id-offset_ );
        if ( par_cg_id != LHAEvent::invalid_id )
          return (unsigned short)ev.at( par_cg_id ).role();
        return findRole( ev, pythia_->event[par_id] );
      }
      return (unsigned short)Particle::UnknownRole;
    }
#endif
  }

  //================================================================================================
  // Custom LHA event definition
  //================================================================================================

#ifdef PYTHIA8
  const double LHAEvent::mp_ = particleproperties::mass( PDG::proton );
  const double LHAEvent::mp2_ = LHAEvent::mp_*LHAEvent::mp_;

  LHAEvent::LHAEvent( const Parameters* params ) :
    LHAup( 3 ), params_( params )
  {
    addProcess( 0, 1., 1., 1.e3 );
    if ( params_ ) {
      setBeamA( (short)params_->kinematics.incoming_beams.first.pdg, params_->kinematics.incoming_beams.first.pz );
      setBeamB( (short)params_->kinematics.incoming_beams.second.pdg, params_->kinematics.incoming_beams.second.pz );
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
  LHAEvent::feedEvent( const Event& ev, bool full, const KinematicsMode& mode )
  {
    const double scale = ev.getOneByRole( Particle::Intermediate ).mass();
    setProcess( 0, 1., scale, constants::alphaEM, constants::alphaQCD );

    const Particle& part1 = ev.getOneByRole( Particle::Parton1 ), &part2 = ev.getOneByRole( Particle::Parton2 );
    const Particle& op1 = ev.getOneByRole( Particle::OutgoingBeam1 ), &op2 = ev.getOneByRole( Particle::OutgoingBeam2 );
    const double q2_1 = -part1.momentum().mass2(), q2_2 = -part2.momentum().mass2();
    const double x1 = q2_1/( q2_1+op1.mass2()-mp2_ ), x2 = q2_2/( q2_2+op2.mass2()-mp2_ );

    unsigned short quark1_id = 0, quark2_id = 0;
    unsigned short quark1_pdgid = part1.integerPdgId(), quark2_pdgid = part2.integerPdgId();

    const Pythia8::Vec4 mom_part1( hadr::momToVec4( part1.momentum() ) ), mom_part2( hadr::momToVec4( part2.momentum() ) );

    if ( !full ) {
      //=============================================================================================
      // incoming partons
      //=============================================================================================

      addCorresp( sizePart(), part1.id() );
      addParticle( quark1_pdgid, -2, quark1_id, 0, 0, 0, mom_part1.px(), mom_part1.py(), mom_part1.pz(), mom_part1.e(), mom_part1.mCalc(), 0., 0. );

      addCorresp( sizePart(), part2.id() );
      addParticle( quark2_pdgid, -2, quark2_id, 0, 0, 0, mom_part2.px(), mom_part2.py(), mom_part2.pz(), mom_part2.e(), mom_part2.mCalc(), 0., 0. );
    }
    else { // full event content (with collinear partons)
      const bool inel1 = ( mode == KinematicsMode::InelasticElastic || mode == KinematicsMode::InelasticInelastic );
      const bool inel2 = ( mode == KinematicsMode::ElasticInelastic || mode == KinematicsMode::InelasticInelastic );

      Pythia8::Vec4 mom_iq1 = mom_part1, mom_iq2 = mom_part2;
      unsigned short colour_index = 501, quark1_colour = 0, quark2_colour = 0;
      //FIXME select quark flavours accordingly
      if ( inel1 ) {
        quark1_pdgid = 2;
        quark1_colour = colour_index++;
        mom_iq1 = hadr::momToVec4( x1*ev.getOneByRole( Particle::IncomingBeam1 ).momentum() );
      }
      if ( inel2 ) {
        quark2_pdgid = 2;
        quark2_colour = colour_index++;
        mom_iq2 = hadr::momToVec4( x2*ev.getOneByRole( Particle::IncomingBeam2 ).momentum() );
      }

      //--- flavour / x value of hard-process initiators
      setIdX( part1.integerPdgId(), part2.integerPdgId(), x1, x2 );

      //===========================================================================================
      // incoming valence quarks
      //===========================================================================================

      quark1_id = sizePart();
      addCorresp( quark1_id, op1.id() );
      addParticle( quark1_pdgid, -1, 0, 0, quark1_colour, 0, mom_iq1.px(), mom_iq1.py(), mom_iq1.pz(), mom_iq1.e(), mom_iq1.mCalc(), 0., 1. );

      quark2_id = sizePart();
      addCorresp( quark2_id, op2.id() );
      addParticle( quark2_pdgid, -1, 0, 0, quark2_colour, 0, mom_iq2.px(), mom_iq2.py(), mom_iq2.pz(), mom_iq2.e(), mom_iq2.mCalc(), 0., 1. );

      //===========================================================================================
      // outgoing valence quarks
      //===========================================================================================

      if ( inel1 ) {
        const Pythia8::Vec4 mom_oq1 = mom_iq1-mom_part1;
        addParticle( quark1_pdgid, 1, quark1_id, quark2_id, quark1_colour, 0, mom_oq1.px(), mom_oq1.py(), mom_oq1.pz(), mom_oq1.e(), mom_oq1.mCalc(), 0., 1. );
      }
      if ( inel2 ) {
        const Pythia8::Vec4 mom_oq2 = mom_iq2-mom_part2;
        addParticle( quark2_pdgid, 1, quark1_id, quark2_id, quark2_colour, 0, mom_oq2.px(), mom_oq2.py(), mom_oq2.pz(), mom_oq2.e(), mom_oq2.mCalc(), 0., 1. );
      }
    }

    //=============================================================================================
    // central system
    //=============================================================================================

    for ( const auto& p : ev.getByRole( Particle::CentralSystem ) ) {
      const auto mothers = p.mothers();
      unsigned short moth1_id = 1, moth2_id = 2;
      if ( !full ) {
        moth1_id = moth2_id = 0;
        if ( mothers.size() > 0 ) {
          const unsigned short moth1_cg_id = *mothers.begin();
          moth1_id = pythiaId( moth1_cg_id );
          if ( moth1_id == invalid_id ) {
            const Particle& moth = ev.at( moth1_cg_id );
            if ( moth.mothers().size() > 0 )
              moth1_id = pythiaId( *moth.mothers().begin() );
            if ( moth.mothers().size() > 1 )
              moth2_id = pythiaId( *moth.mothers().rbegin() );
          }
          if ( mothers.size() > 1 ) {
            const unsigned short moth2_cg_id = *mothers.rbegin();
            moth2_id = pythiaId( moth2_cg_id );
            if ( moth2_id == invalid_id ) {
              const Particle& moth = ev.at( moth2_cg_id );
              moth.dump();
              moth2_id = pythiaId( *moth.mothers().rbegin() );
            }
          }
        }
      }
      const Pythia8::Vec4 mom_part( p.momentum().px(), p.momentum().py(), p.momentum().pz(), p.momentum().energy() );
      addCorresp( sizePart(), p.id() );
      addParticle( p.integerPdgId(), 1, moth1_id, moth2_id, 0, 0, mom_part.px(), mom_part.py(), mom_part.pz(), mom_part.e(), mom_part.mCalc(), 0., 0., 0. );
    }
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
  LHAEvent::cepgenId( unsigned short py_id ) const
  {
    for ( const auto& py_cg : py_cg_corresp_ )
      if ( py_cg.first == py_id )
        return py_cg.second;
    return invalid_id;
  }

  unsigned short
  LHAEvent::pythiaId( unsigned short cg_id ) const
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
    CG_INFO( "LHAEvent:dump" ) << oss.str();
  }
#endif
}
