#include "GenericHadroniser.h"

GenericHadroniser::GenericHadroniser( const std::string& name_ ) :
  fName( name_ ), fHadrons( new std::vector<Particle>() )
{}

GenericHadroniser::~GenericHadroniser()
{
  Debugging( Form("Destructor called" ) );
  
  if ( fHadrons ) delete fHadrons;
}
