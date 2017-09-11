#ifdef PYTHIA8
#include "Pythia8Hadroniser.h"

namespace CepGen
{
  namespace Hadroniser
  {
    Pythia8Hadroniser::Pythia8Hadroniser():
      GenericHadroniser( "Pythia8" ),
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
