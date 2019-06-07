#ifdef PYTHIA8

#include "CepGen/Hadronisers/PythiaEventInterface.h"

#include "CepGen/Core/Exception.h"

#include "CepGen/Parameters.h"

#include "CepGen/Physics/Kinematics.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/Constants.h"

#include "CepGen/Event/Event.h"
#include "CepGen/Event/Particle.h"

namespace Pythia8
{
  /// Convert a CepGen particle momentum into its Pythia8 counterpart
  Vec4
  momToVec4( const cepgen::Particle::Momentum& mom )
  {
    return Vec4( mom.px(), mom.py(), mom.pz(), mom.energy() );
  }

  const double CepGenEvent::mp_ = cepgen::PDG::get()( cepgen::PDG::proton ).mass;
  const double CepGenEvent::mp2_ = CepGenEvent::mp_*CepGenEvent::mp_;

  CepGenEvent::CepGenEvent() :
    LHAup( 3 ),
    inel1_( false ), inel2_( false ), params_( nullptr )
  {}

  void
  CepGenEvent::initialise( const cepgen::Parameters& params )
  {
    params_ = &params;
    const cepgen::KinematicsMode& mode = params_->kinematics.mode;
    inel1_ = ( mode == cepgen::KinematicsMode::InelasticElastic || mode == cepgen::KinematicsMode::InelasticInelastic );
    inel2_ = ( mode == cepgen::KinematicsMode::ElasticInelastic || mode == cepgen::KinematicsMode::InelasticInelastic );

    setBeamA( (short)params_->kinematics.incoming_beams.first.pdg, params_->kinematics.incoming_beams.first.pz );
    setBeamB( (short)params_->kinematics.incoming_beams.second.pdg, params_->kinematics.incoming_beams.second.pz );
    addProcess( 0, params_->integration().result, params_->integration().err_result, 100. );
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
    const double scale = ev[cepgen::Particle::Intermediate][0].mass();
    setProcess( 0, 1., scale, cepgen::constants::ALPHA_EM, cepgen::constants::ALPHA_QCD );

    const auto& part1 = ev[cepgen::Particle::Parton1][0], &part2 = ev[cepgen::Particle::Parton2][0];
    const auto& op1 = ev[cepgen::Particle::OutgoingBeam1][0], &op2 = ev[cepgen::Particle::OutgoingBeam2][0];
    const double q2_1 = -part1.momentum().mass2(), q2_2 = -part2.momentum().mass2();
    const double x1 = q2_1/( q2_1+op1.mass2()-mp2_ ), x2 = q2_2/( q2_2+op2.mass2()-mp2_ );

    unsigned short colour_index = 501;

    const Vec4 mom_part1( momToVec4( part1.momentum() ) ), mom_part2( momToVec4( part2.momentum() ) );

    if ( !full ) {
      //=============================================================================================
      // incoming partons
      //=============================================================================================

      addCorresp( sizePart(), part1.id() );
      addParticle( part1.integerPdgId(), -2, 0, 0, 0, 0, mom_part1.px(), mom_part1.py(), mom_part1.pz(), mom_part1.e(), mom_part1.mCalc(), 0., 0. );

      addCorresp( sizePart(), part2.id() );
      addParticle( part2.integerPdgId(), -2, 0, 0, 0, 0, mom_part2.px(), mom_part2.py(), mom_part2.pz(), mom_part2.e(), mom_part2.mCalc(), 0., 0. );
    }
    else { // full event content (with collinear partons)
      Vec4 mom_iq1 = mom_part1, mom_iq2 = mom_part2;
      unsigned short parton1_id = 0, parton2_id = 0;
      unsigned short parton1_pdgid = part1.integerPdgId(), parton2_pdgid = part2.integerPdgId();
      unsigned short parton1_colour = 0, parton2_colour = 0;
      //FIXME select quark flavours accordingly
      if ( inel1_ ) {
        parton1_pdgid = 2;
        parton1_colour = colour_index++;
        mom_iq1 = momToVec4( x1*ev[cepgen::Particle::IncomingBeam1][0].momentum() );
      }
      if ( inel2_ ) {
        parton2_pdgid = 2;
        parton2_colour = colour_index++;
        mom_iq2 = momToVec4( x2*ev[cepgen::Particle::IncomingBeam2][0].momentum() );
      }

      //--- flavour / x value of hard-process initiators
      setIdX( part1.integerPdgId(), part2.integerPdgId(), x1, x2 );
      setPdf( parton1_pdgid, parton2_pdgid, x1, x2, scale, 0., 0., false );

      //===========================================================================================
      // incoming valence quarks
      //===========================================================================================

      parton1_id = sizePart();
      addCorresp( parton1_id, op1.id() );
      addParticle( parton1_pdgid, -1, 0, 0, parton1_colour, 0, mom_iq1.px(), mom_iq1.py(), mom_iq1.pz(), mom_iq1.e(), mom_iq1.mCalc(), 0., 1. );

      parton2_id = sizePart();
      addCorresp( parton2_id, op2.id() );
      addParticle( parton2_pdgid, -1, 0, 0, parton2_colour, 0, mom_iq2.px(), mom_iq2.py(), mom_iq2.pz(), mom_iq2.e(), mom_iq2.mCalc(), 0., 1. );

      //===========================================================================================
      // outgoing valence quarks
      //===========================================================================================

      if ( inel1_ ) {
        const Vec4 mom_oq1 = mom_iq1-mom_part1;
        addParticle( parton1_pdgid, 1, parton1_id, parton2_id, parton1_colour, 0, mom_oq1.px(), mom_oq1.py(), mom_oq1.pz(), mom_oq1.e(), mom_oq1.mCalc(), 0., 1. );
      }
      if ( inel2_ ) {
        const Vec4 mom_oq2 = mom_iq2-mom_part2;
        addParticle( parton2_pdgid, 1, parton1_id, parton2_id, parton2_colour, 0, mom_oq2.px(), mom_oq2.py(), mom_oq2.pz(), mom_oq2.e(), mom_oq2.mCalc(), 0., 1. );
      }
    }

    //=============================================================================================
    // central system
    //=============================================================================================

    const unsigned short central_colour = colour_index++;
    unsigned short cp_colour = 0, cp_anticolour = 0;
    for ( const auto& p : ev[cepgen::Particle::CentralSystem] ) {
      const auto mothers = p.mothers();
      unsigned short moth1_id = 1, moth2_id = 2;
      if ( !full ) {
        moth1_id = moth2_id = 0;
        if ( !mothers.empty() ) {
          const unsigned short moth1_cg_id = *mothers.begin();
          moth1_id = pythiaId( moth1_cg_id );
          if ( moth1_id == INVALID_ID ) {
            const auto& moth = ev[moth1_cg_id];
            if ( moth.mothers().size() > 0 )
              moth1_id = pythiaId( *moth.mothers().begin() );
            if ( moth.mothers().size() > 1 )
              moth2_id = pythiaId( *moth.mothers().rbegin() );
          }
          if ( mothers.size() > 1 ) {
            const unsigned short moth2_cg_id = *mothers.rbegin();
            moth2_id = pythiaId( moth2_cg_id );
            if ( moth2_id == INVALID_ID ) {
              const auto& moth = ev[moth2_cg_id];
              moth.dump();
              moth2_id = pythiaId( *moth.mothers().rbegin() );
            }
          }
        }
      }
      std::cout << (int)p.pdgId() << std::endl;
      if ( cepgen::PDG::get()( p.pdgId() ).colours > 1 ) {
        if ( p.integerPdgId() > 0 ) //--- particle
          cp_colour = central_colour;
        else //--- anti-particle
          cp_anticolour = central_colour;
      }
      const Vec4 mom_part( momToVec4( p.momentum() ) );
      addCorresp( sizePart(), p.id() );
      addParticle( p.integerPdgId(), 1, moth1_id, moth2_id, cp_colour, cp_anticolour, mom_part.px(), mom_part.py(), mom_part.pz(), mom_part.e(), mom_part.mCalc(), 0., 0., 0. );
    }
  }

  void
  CepGenEvent::setProcess( int id, double xsec, double q2_scale, double alpha_qed, double alpha_qcd )
  {
    LHAup::setProcess( id, xsec, q2_scale, alpha_qed, alpha_qcd );
    py_cg_corresp_.clear();
  }

  unsigned short
  CepGenEvent::cepgenId( unsigned short py_id ) const
  {
    if ( py_cg_corresp_.count( py_id ) == 0 )
      return INVALID_ID;
    return py_cg_corresp_.at( py_id );
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
    py_cg_corresp_[py_id] = cg_id;
  }

  void
  CepGenEvent::dumpCorresp() const
  {
    std::ostringstream oss;
    oss << "List of Pythia ←|→ CepGen particle ids correspondance";
    for ( const auto& py_cg : py_cg_corresp_ )
      oss << "\n\t" << py_cg.first << " <-> " << py_cg.second;
    CG_INFO( "CepGenEvent:dump" ) << oss.str();
  }
}

#endif
