#include <iostream>

#include "CepGen/Generator.h"
#include "CepGen/Cards/ConfigHandler.h"
#include "CepGen/Export/LHEFHandler.h"

#include "HepMC/Version.h"

using namespace std;

/**
 * Main caller for this Monte Carlo generator. Loads the configuration files'
 * variables if set as an argument to this program, else loads a default
 * "LHC-like" configuration, then launches the cross-section computation and
 * the events generation.
 * \author Laurent Forthomme <laurent.forthomme@cern.ch>
 */
int main( int argc, char* argv[] ) {
  CepGen::Generator mg;
  
  if ( argc == 1 ) FatalError( "No config file provided." );

  Debugging( Form( "Reading config file stored in %s", argv[1] ) );
  CepGen::Cards::ConfigHandler card( argv[1] );
  mg.setParameters( card.parameters() );

  // We might want to cross-check visually the validity of our run
  mg.parameters->dump();

  // Let there be cross-section...
  double xsec, err;
  mg.computeXsection( xsec, err );

  CepGen::OutputHandler::LHEFHandler writer( "example.dat" );
  writer.setCrossSection( xsec, err );
  writer.initialise( *mg.parameters );

  // The events generation starts here !
  for ( unsigned int i = 0; i < mg.parameters->generation.maxgen; ++i ) {
    if ( i%1000 == 0 )
      cout << "Generating event #" << i+1 << endl;
    try {
      writer << mg.generateOneEvent();
    } catch ( CepGen::Exception& e ) { e.dump(); }
  }

  return 0;
}

