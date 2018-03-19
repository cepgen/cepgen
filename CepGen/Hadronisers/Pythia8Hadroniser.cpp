#ifdef PYTHIA8
#include "Pythia8Hadroniser.h"

#include "CepGen/Parameters.h"
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
      max_attempts_( params.hadroniser_max_trials )
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
      bool res = pythia_->init();
      if ( !res )
        FatalError( "Failed to initialise the Pythia8 core!\n\t"
                    "See the message above for more details." );
      return res;
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
      //--- start by cleaning up the previous runs leftovers
      pythia_->event.reset();
      py_cg_corresp_.clear();
      cg_py_corresp_.clear();

      //===========================================================================================
      // loop to add the particles
      //===========================================================================================

      for ( unsigned short i = 0; i < ev.numParticles(); ++i ) {
        const Particle& part = ev.getConstById( i );
        if ( part.role() != Particle::CentralSystem )
          continue;

        const Particle::Momentum mom = part.momentum();
        const double mass2 = mom.energy2()-mom.p2(), mass = ( mass2 < 0 ) ? -sqrt( -mass2 ) : sqrt( mass2 );
        Pythia8::Particle py8part( part.integerPdgId(), 23, 0, 0, 0, 0, 0, 0, mom.px(), mom.py(), mom.pz(), mom.energy(), mass );
        unsigned short py_id = pythia_->event.append( py8part );
        // populate the CepGen id <-> Pythia id map
        cg_py_corresp_[part.id()] = py_id;
        py_cg_corresp_[py_id] = part.id();
      }

      //===========================================================================================
      // launch the resonances decays
      //===========================================================================================

      const unsigned short num_py_parts = pythia_->event.size();
      unsigned short num_decays_trials = 0;
      while ( true ) {
        if ( !launchPythia( ev ) )
          return false;
        if ( pythia_->event.size() == num_py_parts )
          return true;

        //--- check if the final state is found
        if ( min_ids_.size() == 0 )
          break;
        std::vector<unsigned short> ids = min_ids_;
        for ( unsigned short i = num_py_parts; i < pythia_->event.size(); ++i ) {
          const Pythia8::Particle& p = pythia_->event[i];
          if ( !p.isFinal() )
            continue;
          std::vector<unsigned short>::iterator it = std::find( ids.begin(), ids.end(), abs( p.id() ) );
          if ( it != ids.end() )
            ids.erase( it );
        }
        if ( ids.size() == 0 )
          break;
        if ( ++num_decays_trials > 5 )
          return false;
      }

      updateEvent( ev, weight );
#endif
      return true;
    }

    bool
    Pythia8Hadroniser::run( Event& ev, double& weight )
    {
      weight = 1.;
#ifndef PYTHIA8
      FatalError( "Pythia8 is not linked to this instance!" );
#else
      const bool full_event = false;
      const double mp = ParticleProperties::mass( Proton ), mp2 = mp*mp;

      //--- start by cleaning up the previous runs leftovers
      pythia_->event.reset();
      py_cg_corresp_.clear();
      cg_py_corresp_.clear();

      //===========================================================================================
      // loop to add the particles
      //===========================================================================================

      unsigned short idx_remn1 = invalid_idx_, idx_remn2 = invalid_idx_;
      for ( unsigned short i = 0; i < ev.numParticles(); ++i ) {
        const Particle& part = ev.getConstById( i );

        const Particle::Momentum mom = part.momentum();
        const double mass2 = mom.energy2()-mom.p2(), mass = ( mass2 < 0 ) ? -sqrt( -mass2 ) : sqrt( mass2 );
        Pythia8::Particle py8part( part.integerPdgId(), 0, 0, 0, 0, 0, 0, 0, mom.px(), mom.py(), mom.pz(), mom.energy(), mass );
        unsigned short py_id = invalid_idx_;
        switch ( part.role() ) {
          case Particle::IncomingBeam1:
          case Particle::IncomingBeam2: {
            if ( !full_event )
              continue;
            py8part.status( -12 );
            py_id = pythia_->event.append( py8part );
          } break;
          case Particle::Parton1:
          case Particle::Parton2:
          case Particle::Parton3: {
            if ( !full_event )
              continue;
            py8part.status( -21 );
            py_id = pythia_->event.append( py8part );
          } break;
          case Particle::Intermediate:
          case Particle::UnknownRole:
            continue;
          case Particle::CentralSystem: {
            if ( !full_event )
              continue;
            py8part.status( -23 ); // outgoing particles of the hardest subprocess
            py_id = pythia_->event.append( py8part );
          } break;
          case Particle::OutgoingBeam1:
          case Particle::OutgoingBeam2: {
            py8part.status(
              ( part.status() == Particle::Unfragmented )
                ? 15
                : 14 // final state proton
            );
//            py8part.scale( 20. );
//            py8part.id( 2212 );
            py_id = pythia_->event.append( py8part );

            if ( part.status() == Particle::Unfragmented ) {
              if ( part.role() == Particle::OutgoingBeam1 )
                idx_remn1 = py_id;
              if ( part.role() == Particle::OutgoingBeam2 )
                idx_remn2 = py_id;
            }
          } break;
        }
        // populate the CepGen id <-> Pythia id map
        if ( py_id != invalid_idx_ ) {
          cg_py_corresp_[part.id()] = py_id;
          py_cg_corresp_[py_id] = part.id();
        }
      }

      //===========================================================================================
      // Particles parentage
      //===========================================================================================

      if ( full_event ) {
        for ( const auto& cg_py : cg_py_corresp_ ) {
          const Particle& part = ev.getConstById( cg_py.first );
          const auto mothers = part.mothers();
          const auto daughters = part.daughters();
          if ( mothers.size() == 0 && daughters.size() == 0 )
            continue;
          { //--- mothers part
            unsigned short id_moth1 = 0, id_moth2 = 0;
            if ( part.role() == Particle::CentralSystem ) {
              const unsigned short id_p1 = ev.getOneByRole( Particle::Parton1 ).id();
              const unsigned short id_p2 = ev.getOneByRole( Particle::Parton2 ).id();
              if ( cg_py_corresp_.count( id_p1 ) > 0 )
                id_moth1 = cg_py_corresp_.at( id_p1 );
              if ( cg_py_corresp_.count( id_p2 ) > 0 )
                id_moth2 = cg_py_corresp_.at( id_p2 );
            }
            else if ( mothers.size() > 0 ) {
              if ( cg_py_corresp_.count( *mothers.begin() ) > 0 )
                id_moth1 = cg_py_corresp_.at( *mothers.begin() );
              if ( mothers.size() > 1 && cg_py_corresp_.count( *mothers.rbegin() ) > 0 )
                id_moth2 = cg_py_corresp_.at( *mothers.rbegin() );
            }
            if ( id_moth2 > id_moth1+1 ) { // id2 <-> id1
              id_moth1 = id_moth1 ^ id_moth2;
              id_moth2 = id_moth1 ^ id_moth2;
              id_moth1 = id_moth1 ^ id_moth2;
            }
            pythia_->event[cg_py.second].mothers( id_moth1, id_moth2 );
          }
          { //--- daughters part
            unsigned short id_daugh1 = 0, id_daugh2 = 0;
            if ( part.role() == Particle::Parton1 || part.role() == Particle::Parton2 ) {
              const auto daugh_gg = ev.getOneByRole( Particle::Intermediate ).daughters();
              if ( daugh_gg.size() == 0 )
                continue;
              if ( cg_py_corresp_.count( *daugh_gg.begin() ) > 0 )
                id_daugh1 = cg_py_corresp_.at( *daugh_gg.begin() );
              if ( daugh_gg.size() > 1 && cg_py_corresp_.count( *daugh_gg.rbegin() ) > 0 )
                id_daugh2 = cg_py_corresp_.at( *daugh_gg.rbegin() );
            }
            else {
              if ( daughters.size() == 0 )
                continue;
              if ( cg_py_corresp_.count( *daughters.begin() ) > 0 )
                id_daugh1 = cg_py_corresp_.at( *daughters.begin() );
              if ( daughters.size() > 1 && cg_py_corresp_.count( *daughters.rbegin() ) > 0 )
                id_daugh2 = cg_py_corresp_.at( *daughters.rbegin() );
            }
            if ( id_daugh2 > id_daugh1+1 ) { // id2 <-> id1
              id_daugh1 = id_daugh1 ^ id_daugh2;
              id_daugh2 = id_daugh1 ^ id_daugh2;
              id_daugh1 = id_daugh1 ^ id_daugh2;
            }
            pythia_->event[cg_py.second].daughters( id_daugh1, id_daugh2 );
            //std::cout << cg_py.second << "|" << abs(pythia_->event[cg_py.second].mCalc() - pythia_->event[cg_py.second].m()) / std::max( 1.0, pythia_->event[cg_py.second].e()) << std::endl;
          }
        }
      }

      //===========================================================================================
      // outgoing remnants massaging
      //===========================================================================================

      if ( idx_remn1 != invalid_idx_ ) {
        const Particle::Momentum& p0 = ev.getOneByRole( Particle::IncomingBeam1 ).momentum();
        const Particle::Momentum& p = ev.getOneByRole( Particle::OutgoingBeam1 ).momentum();
        const double q2 = -( p-p0 ).mass2();
        const double mx2 = p.mass2();
        const double xbj = q2/( q2+mx2-mp2 );
        fragmentState( idx_remn1, xbj );
      }
      if ( idx_remn2 != invalid_idx_ ) {
        const Particle::Momentum& p0 = ev.getOneByRole( Particle::IncomingBeam2 ).momentum();
        const Particle::Momentum& p = ev.getOneByRole( Particle::OutgoingBeam2 ).momentum();
        const double q2 = -( p-p0 ).mass2();
        const double my2 = p.mass2();
        const double xbj = q2/( q2+my2-mp2 );
        fragmentState( idx_remn2, xbj );
      }
      //pythia_->forceHadronLevel();

//pythia_->event.list(true,true);
//return true;

      const unsigned short num_py_parts = pythia_->event.size();
//pythia_->event.list(true,true);
      //ev.dump();
//      exit(0);

      //===========================================================================================
      // launch the hadronisation / resonances decays
      //===========================================================================================

      //std::cout << ":::::" << pythia_->settings.parm("Check:mTolErr") << std::endl;

      if ( !launchPythia( ev ) )
        return false;
      if ( pythia_->event.size() == num_py_parts )
        return true;

//        pythia_->event.list(true,true);
/*      for ( unsigned short i = 0; i < pythia_->event.size(); ++i ) {
        const double err = abs( pythia_->event[i].mCalc()-pythia_->event[i].m() )/std::max( 1.0, pythia_->event[i].e());
        std::cout << ">> " << pythia_->event[i].id() << "::" << err << "::" << (err>pythia_->settings.parm("Check:mTolErr")) << std::endl;
      }*/
//        exit(0);
      //std::cout << "bad" << std::endl;
      //InWarning( "Pythia8 failed to process the event." );
      //pythia_->event.list( true, true );

      // check if something happened in the event processing by Pythia
      // if not, return the event as it is...

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
}

#endif
