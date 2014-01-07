#include "hadroniser.h"

Hadroniser::Hadroniser() :
  _name("undefined")
{
  this->_hadrons = new std::vector<Particle>;
}

Hadroniser::~Hadroniser()
{
  delete [] this->_hadrons;
}

bool
Hadroniser::Hadronise(Particle* part_)
{
  return true;
}
