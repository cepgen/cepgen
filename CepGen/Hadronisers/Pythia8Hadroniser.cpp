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

      /*Pythia8::LHAgenerator gen_info;
      gen_info.name = "CepGen";
      gen_info.version = version();
      pythia_->info.generators->emplace_back( gen_info );*/
#endif
    }

    Pythia8Hadroniser::~Pythia8Hadroniser()
    {}

    void
    Pythia8Hadroniser::setSeed( long long seed )
    {
#ifdef PYTHIA8
      if ( seed == -1ll ) {
        pythia_->readString( "Random:setSeed = off" );
        return;
      }
      pythia_->readString( "Random:setSeed = on" );
      pythia_->readString( Form( "Random:seed = %llu", seed ) );
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
      ids_corresp_.clear();

      //--- first loop to add the particles (and their childs)

      for ( unsigned short i = 0; i < num_before; ++i ) {
        const Particle& part = ev.getConstById( i );
        addParticle( part, ev );
      }

      //--- second loop to specify the particles' parentage

      for ( unsigned short i = 0; i < num_before; ++i ) {
        const Particle& part = ev.getConstById( i );
        const short py_id = ids_corresp_.at( i ); // particle is already there ; retrieve its Pythia8 event id

        Pythia8::Particle& pp = pythia_->event[py_id];
        pp.mother1( ( part.mothers().size() > 0 ) ? ids_corresp_.at( *part.mothers().begin() ) : 0 ); // offset for central system
        pp.mother2( ( part.mothers().size() > 1 ) ? ids_corresp_.at( *part.mothers().rbegin() ) : 0 ); // offset for central system
        pp.daughter1( ( part.daughters().size() > 0 ) ? ids_corresp_.at( *part.daughters().begin() ) : 0 ); // offset for central system
        pp.daughter2( ( part.daughters().size() > 1 ) ? ids_corresp_.at( *part.daughters().rbegin() ) : 0 ); // offset for central system
      }

      ev.dump();
      pythia_->event.list();
      const unsigned short num_py_parts = pythia_->event.size();

      if ( !pythia_->next() ) exit( 0 );//return 0.;

      // check if something happened in the event processing by Pythia
      // if not, return the event as it is...
      if ( pythia_->event.size() == num_py_parts ) return weight;

      pythia_->event.list();

#else
      FatalError( "Pythia8 is not linked to this instance!" );
#endif
      return true;
    }

#ifdef PYTHIA8
    void
    Pythia8Hadroniser::addParticle( const Particle& part, const Event& ev, bool recursive )
    {
      // check if the particle is already inserted
      if ( ids_corresp_.find( part.id() ) != ids_corresp_.end() )
        return;

      const Particle::Momentum mom = part.momentum();
      Pythia8::Particle py8part( part.integerPdgId(), 0, 0, 0, 0, 0, 0, 0, mom.px(), mom.py(), mom.pz(), mom.energy(), part.mass() );
      switch ( part.role() ) {
        case Particle::IncomingBeam1:
        case Particle::IncomingBeam2: {
          py8part.status( -12 );
        } break;
        case Particle::Parton1:
        case Particle::Parton2:
        case Particle::Parton3: {
          //py8part.status( -13 );
        } break;
        case Particle::Intermediate: {
          py8part.id( ev.getOneByRole( Particle::Parton1 ).integerPdgId() );
          //py8part.status( -14 );
        }
        case Particle::CentralSystem: {
          if ( pythia_->particleData.canDecay( py8part.id() ) ) {
            py8part.status( 93 );
            py8part.tau( pythia_->particleData.tau0( py8part.id() ) );
            //part.setStatus( Particle::Resonance );
          }
        } break;
        case Particle::OutgoingBeam1:
        case Particle::OutgoingBeam2: {
          if ( py8part.id() == 2 ) { // dissociative proton
            py8part.id( 9902210 );
            //py8part.status( 15 );
            py8part.status( 63 );
          }
        } break;
        default: break;
      }
      ids_corresp_.insert( std::pair<short,short>( part.id(), pythia_->event.append( py8part ) ) );

      if ( recursive ) {
        const ParticlesIds daughters = part.daughters();
        for ( ParticlesIds::const_iterator it = daughters.begin(); it != daughters.end(); ++it ) {
          addParticle( ev.getConstById( *it ), ev, false );
        }
      }
    }

    double
    Pythia8Hadroniser::decay( const Particle& part, Event& ev )
    {
      double br_ratio = 0.;

      // check if the particle can be decayed
      if ( part.status() != Particle::Undecayed ) return 0.;
      if ( !pythia_->particleData.canDecay( part.pdgId() ) ) return 0.;

      // check if the particle has a valid kinematics
      const Particle::Momentum mom = part.momentum();
      if ( mom == Particle::Momentum() ) return 0.; // FIXME probably during warm-up...

      pythia_->event.reset();
      Pythia8::Particle py8part( part.integerPdgId(), 93, 0, 0, 0, 0, 0, 0, mom.px(), mom.py(), mom.pz(), mom.energy(), part.mass() );
      py8part.tau( pythia_->particleData.tau0( part.pdgId() ) );
      pythia_->event.append( py8part );

      const unsigned short num_before = pythia_->event.size(); // number of particles before the decay
      if ( !pythia_->next() ) return 0.; // launch the decay

      // map { "pythia id" -> "cepgen id" }
      std::map<unsigned short, unsigned short> py_cg_corresp = { { 1, part.id() } };

      const unsigned short num_after = pythia_->event.size(); // number of particles after the decay
      if ( num_before == num_after ) return 0.; // nothing happened... so no decay!

      // check the branching fraction for the given decay
      br_ratio = pythia_->event[1].particleDataEntry().pickChannel().bRatio();
      if ( br_ratio == 0. ) return 0.;

      for ( unsigned short i = 2; i < pythia_->event.size(); ++i ) { // skip the initial system
        const Pythia8::Particle p = pythia_->event[i];
        const std::vector<int> mothers = p.motherList();
        if ( mothers.size() == 0 ) continue;
        Particle& op = ev.addParticle( Particle::CentralSystem );
        py_cg_corresp[i] = op.id();

        op.setPdgId( ( ParticleCode )abs( p.id() ), p.charge() );
        if ( p.isFinal() ) op.setStatus( Particle::FinalState );
        else op.setStatus( Particle::Propagator );

        op.setMomentum( Particle::Momentum( p.px(), p.py(), p.pz(), p.e() ) );
        for ( std::vector<int>::const_iterator moth = mothers.begin(); moth != mothers.end(); ++moth ) {
          if ( py_cg_corresp.count( *moth ) == 0 ) {
            FatalError( Form( "Particle with id=%d was not found in the event content!", *moth ) );
          }
          op.addMother( ev.getById( py_cg_corresp[*moth] ) );
        }
      }
      return br_ratio;
    }
#endif

  }
}

#endif
