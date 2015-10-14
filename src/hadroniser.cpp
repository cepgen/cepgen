#include "hadroniser.h"

Hadroniser::Hadroniser() :
  _name("undefined")
{
  _hadrons = new std::vector<Particle>();
}

Hadroniser::~Hadroniser()
{
  delete _hadrons;
#ifdef DEBUG
  PrintDebug(Form("Destructor called"));
#endif
}

