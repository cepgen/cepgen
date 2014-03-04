#include "pythia8hadroniser.h"

Pythia8Hadroniser::Pythia8Hadroniser()
{
}

Pythia8Hadroniser::~Pythia8Hadroniser()
{
}

bool
Pythia8Hadroniser::Hadronise(Event* ev_)
{
  ev_->Dump();

  //Pythia8::Pythia py;
  //Pythia8::StringFragmentation string;
  //Pythia8::BeamParticle bp;

  Particles::iterator p;
  Particles part;

  part = ev_->GetParticles();

  for (p=part.begin(); p!=part.end(); p++) {
    std::cout << "--> " << (*p)->pdgId << ", " << (*p)->status << std::endl;
    if ((*p)->status==3) {
      
    }
  }

  //Pythia8::ColSinglet cs();
  Pythia8::Pythia py;
  py.readString("");
  

  return true;
}
