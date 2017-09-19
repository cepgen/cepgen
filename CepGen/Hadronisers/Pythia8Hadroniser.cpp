#ifdef PYTHIA8
#include "Pythia8Hadroniser.h"

namespace CepGen
{
  namespace Hadroniser
  {
    Pythia8Hadroniser::Pythia8Hadroniser():
      GenericHadroniser( "pythia8" ),
      pythia_( std::unique_ptr<Pythia8::Pythia>( new Pythia8::Pythia ) )
    {}

    Pythia8Hadroniser::~Pythia8Hadroniser()
    {}

    void
    Pythia8Hadroniser::setSeed( long long seed )
    {
      if ( seed == -1ll ) {
        pythia_->readString( "Random:setSeed = off" );
        return;
      }
      pythia_->readString( "Random:setSeed = on" );
      pythia_->readString( Form( "Random:seed = %llu", seed ) );
    }

    bool
    Pythia8Hadroniser::hadronise( const Particle* part )
    {
      return true;
    }

    bool
    Pythia8Hadroniser::hadronise( Event* ev )
    {
      //--- first start by cleaning up the previous runs leftovers

      pythia_->event.reset();
      std::cout << ">>>>> before hadronisation" << std::endl;
      ev->dump();

      //--- outgoing proton remnants fragmentation

      Particles& op1 = ev->getByRole( Particle::OutgoingBeam1 ), &op2 = ev->getByRole( Particle::OutgoingBeam2 );
      for ( Particles::const_iterator p_it = op1.begin(); p_it != op1.end(); ++p_it ) {
        //... excited proton fragmentation
      }
      for ( Particles::const_iterator p_it = op2.begin(); p_it != op2.end(); ++p_it ) {
        //... excited proton fragmentation
      }

      //--- central system hadronisation/decay/...

      const ParticlesIds cs = ev->getIdsByRole( Particle::CentralSystem );
      for ( ParticlesIds::iterator p_it = cs.begin(); p_it != cs.end(); ++p_it ) {
        const Particle& part = ev->getById( *p_it );

        // check if the particle can be decayed
        if ( part.status() != Particle::Undecayed ) continue;
        if ( !pythia_->particleData.canDecay( part.pdgId() ) ) continue;

        // check if the particle has a valid kinematics
        const Particle::Momentum mom = part.momentum();
        if ( mom == Particle::Momentum() ) continue; // FIXME probably during warm-up...

        pythia_->event.reset();
        Pythia8::Particle py8part( part.integerPdgId(), 93, 0, 0, 0, 0, 0, 0, mom.px(), mom.py(), mom.pz(), mom.energy(), part.mass() );
        py8part.tau( pythia_->particleData.tau0( part.pdgId() ) );
        pythia_->event.append( py8part );

        const unsigned short num_before = pythia_->event.size(); // number of particles before the decay
        pythia_->next(); // launch the decay

        // map { "pythia id" -> "cepgen id" }
        std::map<unsigned short, unsigned short> py_cg_corresp = { { 1, part.id() } };

        const unsigned short num_after = pythia_->event.size(); // number of particles after the decay
        if ( num_before == num_after ) continue; // nothing happened... so no decay!

        for ( unsigned short i = 2; i < pythia_->event.size(); ++i ) { // skip the initial system
          const Pythia8::Particle p = pythia_->event[i];
          const std::vector<int> mothers = p.motherList();
          if ( mothers.size() == 0 ) continue;
          Particle& op = ev->addParticle( Particle::CentralSystem );
          py_cg_corresp[i] = op.id();

//std::cout << "dumping ids correspondence (" << *p_it << ")" << std::endl; for ( const auto& p : py_cg_corresp ) { std::cout << p.first << "--->" << p.second << std::endl; } std::cout << "-------enddump-------" << std::endl;
          op.setPdgId( ( Particle::ParticleCode )abs( p.id() ), p.charge() );
          if ( p.isFinal() ) op.setStatus( Particle::FinalState );
          else op.setStatus( Particle::Propagator );

          op.setMomentum( Particle::Momentum( p.px(), p.py(), p.pz(), p.e() ) );
          for ( std::vector<int>::const_iterator moth = mothers.begin(); moth != mothers.end(); ++moth ) {
            if ( py_cg_corresp.count( *moth ) == 0 ) { FatalError( Form( "Particle with id=%d was not found in the event content!", *moth ) ); }
            op.addMother( ev->getById( py_cg_corresp[*moth] ) );
          }
        }
        ev->getById( *p_it ).setStatus( Particle::Resonance );
        //std::cout << part.id() << "->" << part.integerPdgId() << "->" << part.status() << std::endl;
      }

      ev->dump();

      return true;
    }

    bool
    Pythia8Hadroniser::prepareHadronisation( Event *ev )
    {
      init();
      return true;
    }
  }
}

#endif
