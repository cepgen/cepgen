#include "GenericHadroniser.h"

using namespace CepGen::Hadroniser;

GenericHadroniser::GenericHadroniser( const char* name ) :
  name_( name ), hadrons_( 0 )
{}

GenericHadroniser::~GenericHadroniser()
{
  Debugging( Form("Destructor called" ) );
  
  if ( hadrons_ ) delete hadrons_;
}
