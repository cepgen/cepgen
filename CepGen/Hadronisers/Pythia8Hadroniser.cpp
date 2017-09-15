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

    bool
    Pythia8Hadroniser::hadronise( const Particle* part )
    {
      return true;
    }

    bool
    Pythia8Hadroniser::hadronise( Event* ev )
    {
      //--- outgoing proton remnants fragmentation

      Particles& op1 = ev->getByRole( Particle::OutgoingBeam1 ), &op2 = ev->getByRole( Particle::OutgoingBeam2 );
      for ( Particles::const_iterator p_it = op1.begin(); p_it != op1.end(); ++p_it ) {
        //... excited proton fragmentation
      }
      for ( Particles::const_iterator p_it = op2.begin(); p_it != op2.end(); ++p_it ) {
        //... excited proton fragmentation
      }

      //--- central system hadronisation/decay/...

      Particles& cs = ev->getByRole( Particle::CentralSystem );
      for ( Particles::const_iterator p_it = cs.begin(); p_it != cs.end(); ++p_it ) {
        if ( p_it->status() != Particle::Undecayed ) continue;
        p_it->dump();
      }

      return true;
    }

    bool
    Pythia8Hadroniser::prepareHadronisation( Event *ev )
    {
      pythia_->init();
      return true;
    }
  }
}

#endif
