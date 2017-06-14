#ifdef PYTHIA8
#include "Pythia8Hadroniser.h"

Pythia8Hadroniser::Pythia8Hadroniser():
  GenericHadroniser( "Pythia8" )
{
  fPy = new Pythia8::Pythia;
}

Pythia8Hadroniser::~Pythia8Hadroniser()
{
  delete fPy;
}

bool
Pythia8Hadroniser::Hadronise( const Particle* part )
{
  return true;
}

bool
Pythia8Hadroniser::Hadronise( Event* ev )
{
  return true;
}

bool
Pythia8Hadroniser::PrepareHadronisation( Event *ev )
{
  fPy->init();
  return true;
}

#endif
