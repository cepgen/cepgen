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
      GenericHadroniser( "pythia8" ), max_attempts_( params.hadroniser_max_trials )
    {
#ifdef PYTHIA8
      pythia_.reset( new Pythia8::Pythia );
      std::cout << params.kinematics.inpdg.first << "\t" << params.kinematics.inpdg.second << std::endl;
      readString( Form( "Beams:idA = %d", params.kinematics.inpdg.first ) );
      readString( Form( "Beams:idB = %d", params.kinematics.inpdg.second ) );
      readString( Form( "Beams:eCM = %.2f", params.kinematics.sqrtS() ) );
#endif
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
    Pythia8Hadroniser::hadronise( Event& ev, double& weight, bool proton_fragment )
    {
      weight = 1.;
#ifdef PYTHIA8
      //--- start by cleaning up the previous runs leftovers
      const double mp = ParticleProperties::mass( Proton ), mp2 = mp*mp;

      pythia_->event.reset();
      /*const double sqrt_s = ev.cmEnergy();
      pythia_->event[0].e( sqrt_s );
      pythia_->event[0].m( sqrt_s );*/
      Particle::Momentum p_out = ev.getOneByRole( Particle::OutgoingBeam1 ).momentum()
                               + ev.getOneByRole( Particle::OutgoingBeam2 ).momentum();
      for ( const auto& p : ev.getByRole( Particle::CentralSystem ) )
        p_out += p.momentum();
      //std::cout << p_out.pz() << "|" << p_out.mass() << "|" << p_out.energy() << std::endl;
      pythia_->event[0].e( p_out.energy() );
      pythia_->event[0].px( p_out.px() );
      pythia_->event[0].py( p_out.py() );
      pythia_->event[0].pz( p_out.pz() );
      pythia_->event[0].m( p_out.mass() );

      const unsigned short num_before = ev.numParticles();
      std::map<short,short> py_cg_corresp, cg_py_corresp;

      //--- loop to add the particles

      unsigned short idx_remn1 = invalid_idx_, idx_remn2 = invalid_idx_;
      for ( unsigned short i = 0; i < num_before; ++i ) {
        const Particle& part = ev.getConstById( i );

        const Particle::Momentum mom = part.momentum();
        Pythia8::Particle py8part( part.integerPdgId(), 0, 0, 0, 0, 0, 0, 0, mom.px(), mom.py(), mom.pz(), mom.energy(), part.mass() );
        unsigned short py_id = invalid_idx_;
        switch ( part.role() ) {
          case Particle::IncomingBeam1:
          case Particle::IncomingBeam2:
          case Particle::Parton1:
          case Particle::Parton2:
          case Particle::Parton3:
          case Particle::Intermediate:
          case Particle::UnknownRole:
            continue;
          case Particle::CentralSystem: {
            if ( pythia_->particleData.canDecay( py8part.id() )
              || pythia_->particleData.mayDecay( py8part.id() ) )
              py8part.status( 93 );
            else
              py8part.status( 23 ); // outgoing particles of the hardest subprocess
            py_id = pythia_->event.append( py8part );
          } break;
          case Particle::OutgoingBeam1:
          case Particle::OutgoingBeam2: {
            py8part.id( 2212 );
            py8part.status( ( proton_fragment && part.status() == Particle::Unfragmented )
              ? 15
              : 14 // final state proton
            );
            py_id = pythia_->event.append( py8part );

            if ( proton_fragment && part.status() == Particle::Unfragmented ) {
              if ( part.role() == Particle::OutgoingBeam1 ) idx_remn1 = py_id;
              if ( part.role() == Particle::OutgoingBeam2 ) idx_remn2 = py_id;
            }
          } break;
        }
        // populate the CepGen id <-> Pythia id map
        if ( py_id != invalid_idx_ ) {
          cg_py_corresp[part.id()] = py_id;
          py_cg_corresp[py_id] = part.id();
        }
      }
      if ( proton_fragment ) {
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
      }

      const unsigned short num_py_parts = pythia_->event.size();
//pythia_->event.list(true,true);
      //ev.dump();
//      exit(0);

      //std::cout << ":::::" << pythia_->settings.parm("Check:mTolErr") << std::endl;

      bool success = false;
      ev.num_hadronisation_trials = 0;
      while ( !success ) {
        success = pythia_->next();
        if ( proton_fragment ){
          std::cout << success << "attempt " << ev.num_hadronisation_trials << " / " << max_attempts_ << std::endl;
          pythia_->event.list(true,true);
      }
        if ( ++ev.num_hadronisation_trials > max_attempts_ )
          return false;
      }
      /*for ( unsigned short i = 0; i < pythia_->event.size(); ++i ) {
        const double err = abs( pythia_->event[i].mCalc()-pythia_->event[i].m() )/std::max( 1.0, pythia_->event[i].e());
        std::cout << ">> " << pythia_->event[i].id() << "::" << err << "::" << (err>pythia_->settings.parm("Check:mTolErr")) << std::endl;
      }*/
      //std::cout << "bad" << std::endl;
      //InWarning( "Pythia8 failed to process the event." );
      //pythia_->event.list( true, true );

      // check if something happened in the event processing by Pythia
      // if not, return the event as it is...
      if ( pythia_->event.size() == num_py_parts )
        return true;

      for ( unsigned short i = 1; i < pythia_->event.size(); ++i ) { // skip the central system
        const Pythia8::Particle& p = pythia_->event[i];
        std::map<short,short>::const_iterator it = py_cg_corresp.find( i );
        if ( it != py_cg_corresp.end() ) { // the particle is already in the event content
          if ( p.daughterList().size() > 0 && i != idx_remn1 && i != idx_remn2 ) {
            weight *= p.particleDataEntry().pickChannel().bRatio();
            ev.getById( it->second ).setStatus( Particle::Resonance );
          }
        }
        else { // the particle was not yet included in the CepGen event
          const std::vector<int> mothers = p.motherList();
          if ( mothers.size() == 0 ) continue; // isolated particle

          Particle::Role role = Particle::CentralSystem;
          std::map<short,short>::const_iterator it_m = py_cg_corresp.find( mothers[0] );
          if ( it_m != py_cg_corresp.end() ) {
            const Particle& moth = ev.getById( it_m->second );
            if ( mothers[0] == idx_remn1 || moth.role() == Particle::OutgoingBeam1 )
              role = Particle::OutgoingBeam1;
            else if ( mothers[0] == idx_remn2 || moth.role() == Particle::OutgoingBeam2 )
              role = Particle::OutgoingBeam2;
          }

          Particle& op = ev.addParticle( role );
          cg_py_corresp[op.id()] = i;
          py_cg_corresp[i] = op.id();

          op.setPdgId( static_cast<ParticleCode>( abs( p.id() ) ), p.charge() );
          op.setStatus( p.isFinal()
            ? Particle::FinalState
            : Particle::Propagator
          );
          op.setMomentum( Particle::Momentum( p.px(), p.py(), p.pz(), p.e() ) );
          for ( std::vector<int>::const_iterator moth = mothers.begin(); moth != mothers.end(); ++moth ) {
            if ( *moth != 0 && py_cg_corresp.count( *moth ) == 0 )
              FatalError( Form( "Particle with id=%d was not found in the event content!", *moth ) );
            op.addMother( ev.getById( py_cg_corresp[*moth] ) );
          }
        }
      }
#else
      FatalError( "Pythia8 is not linked to this instance!" );
#endif
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
  }
}

#endif
