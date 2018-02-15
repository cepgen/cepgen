#ifdef PYTHIA8
#include "Pythia8Hadroniser.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Event/Particle.h"
#include "CepGen/Version.h"

namespace CepGen
{
  namespace Hadroniser
  {
    Pythia8Hadroniser::Pythia8Hadroniser() :
      GenericHadroniser( "pythia8" )
    {
#ifdef PYTHIA8
      pythia_.reset( new Pythia8::Pythia );

      //--- start by disabling some unnecessary output
      pythia_->readString( "Next:numberCount = 0" );
#endif
    }

    Pythia8Hadroniser::~Pythia8Hadroniser()
    {}

    bool
    Pythia8Hadroniser::init()
    {
      bool res = pythia_->init();
      if ( !res ) {
        FatalError( "Failed to initialise the Pythia8 core!\n\t"
                    "See the message above for more details." );
      }
      return res;
    }

    void
    Pythia8Hadroniser::readString( const char* param )
    {
      if ( !pythia_->readString( param ) ) {
        FatalError( Form( "The Pythia8 core failed to parse the following setting:\n\t%s", param ) );
      }
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
    Pythia8Hadroniser::hadronise( Event& ev, double& weight )
    {
      weight = 1.;
#ifdef PYTHIA8
      //--- start by cleaning up the previous runs leftovers

      const double sqrt_s = ev.cmEnergy();

      pythia_->event.reset();
      pythia_->event[0].e( sqrt_s );
      pythia_->event[0].m( sqrt_s );

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
            if ( !pythia_->particleData.canDecay( py8part.id() ) )
              continue;
            if ( !pythia_->particleData.mayDecay( py8part.id() ) )
              continue;
            py8part.status( 93 );
            py_id = pythia_->event.append( py8part );
          } break;
          case Particle::OutgoingBeam1:
          case Particle::OutgoingBeam2: {
            py8part.id( 2212 );
            py8part.status( ( part.status() == Particle::Unfragmented )
              ? 93
              : 23 // final state proton
            );
            py_id = pythia_->event.append( py8part );

            if ( part.status() == Particle::Unfragmented ) {
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
      if ( idx_remn1 != invalid_idx_ ) fragmentState( idx_remn1 );
      if ( idx_remn2 != invalid_idx_ ) fragmentState( idx_remn2 );

      const unsigned short num_py_parts = pythia_->event.size();

      if ( !pythia_->next() ) {
        InWarning( "Pythia8 failed to process the event." );
        return false;
      }

//pythia_->event.list(true,true);
//exit(0);

      // check if something happened in the event processing by Pythia
      // if not, return the event as it is...
      if ( pythia_->event.size() == num_py_parts ) return true;

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
      //ev.dump();
#else
      FatalError( "Pythia8 is not linked to this instance!" );
#endif
      return true;
    }

    void
    Pythia8Hadroniser::fragmentState( unsigned short idx )
    {
      const Pythia8::Particle& remn = pythia_->event[idx];
      const Pythia8::Particle& system = pythia_->event[remn.mother1()];
      const short sign = remn.pz()/fabs( remn.pz() );
      const double px_dq = 0., py_dq = 0., pz_dq = sign * system.e()/2.;
      const double rnd = 1./RAND_MAX * rand();
      unsigned short pdg_q = 0, pdg_dq = 0;
      if      ( rnd < 1./9. ) { pdg_q = 1; pdg_dq = 2203; }
      else if ( rnd < 5./9. ) { pdg_q = 2; pdg_dq = 2101; }
      else                    { pdg_q = 2; pdg_dq = 2103; }
      const double m_q = pythia_->particleData.m0( pdg_q ), m_dq = pythia_->particleData.m0( pdg_dq );
      Pythia8::Particle diquark( pdg_dq, 21, idx, 0, 0, 0, 0, 100+idx, px_dq, py_dq, pz_dq, 0., m_dq );
      Pythia8::Particle quark( pdg_q, 21, idx, 0, 0, 0, 100+idx, 0, remn.px()-px_dq, remn.py()-py_dq, remn.pz()-pz_dq, 0., m_q );
      diquark.e( diquark.eCalc() );
      quark.e( quark.eCalc() );
      const unsigned short id_dq = pythia_->event.append( diquark ), id_q = pythia_->event.append( quark );
      pythia_->event[idx].daughter1( id_dq );
      pythia_->event[idx].daughter2( id_q );
      //pythia_->event[idx].status( -pythia_->event[idx].status() );
      pythia_->event[idx].status( -15 );
    }
  }
}

#endif
