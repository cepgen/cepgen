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
  std::cout << __PRETTY_FUNCTION__ << " [DEBUG] Destructor called" << std::endl;
#endif
}

