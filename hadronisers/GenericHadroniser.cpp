#include "GenericHadroniser.h"

GenericHadroniser::GenericHadroniser( const char* name_ ) :
  fName( name_ ), fHadrons( new Particles() )
{}

GenericHadroniser::~GenericHadroniser()
{
  Debugging( Form("Destructor called" ) );
  
  if ( fHadrons ) delete fHadrons;
}
