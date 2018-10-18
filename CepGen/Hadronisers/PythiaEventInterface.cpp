#include "CepGen/Hadronisers/PythiaEventInterface.h"

#include "CepGen/Core/Exception.h"

#include "CepGen/Parameters.h"

#include "CepGen/Physics/Kinematics.h"
#include "CepGen/Physics/ParticleProperties.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/Constants.h"

#include "CepGen/Event/Event.h"
#include "CepGen/Event/Particle.h"

#ifdef PYTHIA8
namespace Pythia8
{
  /// Convert a CepGen particle momentum into its Pythia8 counterpart
  Vec4
  momToVec4( const cepgen::Particle::Momentum& mom )
  {
    return Vec4( mom.px(), mom.py(), mom.pz(), mom.energy() );
  }

  const double CepGenEvent::mp_ = cepgen::particleproperties::mass( cepgen::PDG::proton );
  const double CepGenEvent::mp2_ = CepGenEvent::mp_*CepGenEvent::mp_;

  CepGenEvent::CepGenEvent() :
    LHAup( 3 )
  {}

  void
  CepGenEvent::initialise( const cepgen::Parameters& params )
  {
    setBeamA( (short)params.kinematics.incoming_beams.first.pdg, params.kinematics.incoming_beams.first.pz );
    setBeamB( (short)params.kinematics.incoming_beams.second.pdg, params.kinematics.incoming_beams.second.pz );
    addProcess( 0, params.integrator.result, params.integrator.err_result, 100. );
  }

  void
  CepGenEvent::setCrossSection( int id, double xsec, double xsec_err )
  {
    setXSec( id, xsec );
    setXErr( id, xsec_err );
    //listInit();
  }

  void
  CepGenEvent::feedEvent( const cepgen::Event& ev, bool full )
  {
    const double scale = ev.getOneByRole( cepgen::Particle::Intermediate ).mass();
    setProcess( 0, 1., scale, cepgen::constants::ALPHA_EM, cepgen::constants::ALPHA_QCD );

    const auto& part1 = ev.getOneByRole( cepgen::Particle::Parton1 ), &part2 = ev.getOneByRole( cepgen::Particle::Parton2 );
    const auto& op1 = ev.getOneByRole( cepgen::Particle::OutgoingBeam1 ), &op2 = ev.getOneByRole( cepgen::Particle::OutgoingBeam2 );
    const double q2_1 = -part1.momentum().mass2(), q2_2 = -part2.momentum().mass2();
    const double x1 = q2_1/( q2_1+op1.mass2()-mp2_ ), x2 = q2_2/( q2_2+op2.mass2()-mp2_ );

    unsigned short quark1_id = 0, quark2_id = 0;
    unsigned short quark1_pdgid = part1.integerPdgId(), quark2_pdgid = part2.integerPdgId();

    const Pythia8::Vec4 mom_part1( Pythia8::momToVec4( part1.momentum() ) ), mom_part2( Pythia8::momToVec4( part2.momentum() ) );

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
      const cepgen::KinematicsMode& mode = params_->kinematics.mode;
      const bool inel1 = ( mode == cepgen::KinematicsMode::InelasticElastic || mode == cepgen::KinematicsMode::InelasticInelastic );
      const bool inel2 = ( mode == cepgen::KinematicsMode::ElasticInelastic || mode == cepgen::KinematicsMode::InelasticInelastic );

      Pythia8::Vec4 mom_iq1 = mom_part1, mom_iq2 = mom_part2;
      unsigned short colour_index = 501, quark1_colour = 0, quark2_colour = 0;
      //FIXME select quark flavours accordingly
      if ( inel1 ) {
        quark1_pdgid = 2;
        quark1_colour = colour_index++;
        mom_iq1 = Pythia8::momToVec4( x1*ev.getOneByRole( cepgen::Particle::IncomingBeam1 ).momentum() );
      }
      if ( inel2 ) {
        quark2_pdgid = 2;
        quark2_colour = colour_index++;
        mom_iq2 = Pythia8::momToVec4( x2*ev.getOneByRole( cepgen::Particle::IncomingBeam2 ).momentum() );
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

    for ( const auto& p : ev[cepgen::Particle::CentralSystem] ) {
      const auto mothers = p.mothers();
      unsigned short moth1_id = 1, moth2_id = 2;
      if ( !full ) {
        moth1_id = moth2_id = 0;
        if ( mothers.size() > 0 ) {
          const unsigned short moth1_cg_id = *mothers.begin();
          moth1_id = pythiaId( moth1_cg_id );
          if ( moth1_id == INVALID_ID ) {
            const auto& moth = ev.at( moth1_cg_id );
            if ( moth.mothers().size() > 0 )
              moth1_id = pythiaId( *moth.mothers().begin() );
            if ( moth.mothers().size() > 1 )
              moth2_id = pythiaId( *moth.mothers().rbegin() );
          }
          if ( mothers.size() > 1 ) {
            const unsigned short moth2_cg_id = *mothers.rbegin();
            moth2_id = pythiaId( moth2_cg_id );
            if ( moth2_id == INVALID_ID ) {
              const auto& moth = ev.at( moth2_cg_id );
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

  void
  CepGenEvent::setProcess( int id, double xsec, double q2_scale, double alpha_qed, double alpha_qcd )
  {
    LHAup::setProcess( id, xsec, q2_scale, alpha_qed, alpha_qcd );
    py_cg_corresp_.clear();
  }

  void
  CepGenEvent::addComments( const std::string& comments )
  {
    osLHEF << comments;
  }

  unsigned short
  CepGenEvent::cepgenId( unsigned short py_id ) const
  {
    for ( const auto& py_cg : py_cg_corresp_ )
      if ( py_cg.first == py_id )
        return py_cg.second;
    return INVALID_ID;
  }

  unsigned short
  CepGenEvent::pythiaId( unsigned short cg_id ) const
  {
    for ( const auto& py_cg : py_cg_corresp_ )
      if ( py_cg.second == cg_id )
        return py_cg.first;
    return INVALID_ID;
  }

  void
  CepGenEvent::addCorresp( unsigned short py_id, unsigned short cg_id )
  {
    py_cg_corresp_.emplace_back( py_id, cg_id );
  }

  void
  CepGenEvent::dumpCorresp() const
  {
    std::ostringstream oss;
    oss << "List of Pythia <-> CepGen particle ids correspondance";
    for ( const auto& py_cg : py_cg_corresp_ )
      oss << "\n\t" << py_cg.first << " <-> " << py_cg.second;
    CG_INFO( "CepGenEvent:dump" ) << oss.str();
  }
}
#endif

