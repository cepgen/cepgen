#include "GenericHadroniser.h"

GenericHadroniser::GenericHadroniser(std::string name_) :
  fName(name_), fHadrons(new std::vector<Particle>())
{}

GenericHadroniser::~GenericHadroniser()
{
  Debug(Form("Destructor called"));
  
  delete fHadrons;
}

