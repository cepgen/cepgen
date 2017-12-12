#ifdef PYTHIA8
#include "Pythia8Hadroniser.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Event/Particle.h"

namespace CepGen
{
  namespace Hadroniser
  {
    Pythia8Hadroniser::Pythia8Hadroniser() :
      GenericHadroniser( "pythia8" )
    {
#ifdef PYTHIA8
      pythia_ = std::unique_ptr<Pythia8::Pythia>( new Pythia8::Pythia );

      //--- start by disabling some unnecessary output
      pythia_->readString( "Next:numberCount = 0" );
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
      //--- first start by cleaning up the previous runs leftovers

      pythia_->event.reset();

      //--- then retrieve the collections
      const ParticlesIds cs = ev.getIdsByRole( Particle::CentralSystem ), op1 = ev.getIdsByRole( Particle::OutgoingBeam1 ), op2 = ev.getIdsByRole( Particle::OutgoingBeam2 );

      //--- outgoing proton remnants fragmentation

      for ( ParticlesIds::iterator p_it = op1.begin(); p_it != op1.end(); ++p_it ) {
        // excited proton fragmentation
        const Particle& part = ev.getById( *p_it );
        if ( !fragment( part, ev ) ) continue;
        ev.getById( *p_it ).setStatus( Particle::Undefined );
      }
      for ( ParticlesIds::iterator p_it = op2.begin(); p_it != op2.end(); ++p_it ) {
        // excited proton fragmentation
        const Particle& part = ev.getById( *p_it );
        if ( !fragment( part, ev ) ) continue;
        ev.getById( *p_it ).setStatus( Particle::Undefined );
      }

      //--- central system hadronisation/decay/...

      for ( ParticlesIds::iterator p_it = cs.begin(); p_it != cs.end(); ++p_it ) {
        const Particle& part = ev.getById( *p_it );
        const double br = decay( part, ev );
        if ( br == 0. ) continue;
        weight *= br;
        ev.getById( *p_it ).setStatus( Particle::Resonance );
      }
#else
      FatalError( "Pythia8 is not linked to this instance!" );
#endif
      return true;
    }

#ifdef PYTHIA8
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

    bool
    Pythia8Hadroniser::fragment( const Particle& part, Event& ev )
    {
      // check if the particle can be decayed
      if ( part.status() != Particle::Unfragmented ) return false;
      //if ( !pythia_->particleData.canDecay( part.pdgId() ) ) return false;

      // check if the particle has a valid kinematics
      const Particle::Momentum mom = part.momentum();
      if ( mom == Particle::Momentum() ) return false; // FIXME probably during warm-up...

      pythia_->event.reset();


      return true;
    }
#endif

  }
}

#endif
