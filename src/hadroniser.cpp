#include "hadroniser.h"

Hadroniser::Hadroniser(std::string name_) :
  fName(name_), fHadrons(new std::vector<Particle>())
{}

Hadroniser::~Hadroniser()
{
  Debug(Form("Destructor called"));
  
  delete fHadrons;
}

