#include "CepGen/Cards/PythonHandler.h"
#include "CepGen/Generator.h"
#include "CepGen/Export/HepMCHandler.h"
#include "CepGen/Core/utils.h"

#include "HepMC/Version.h"

#include <iostream>

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

  if ( argc == 1 )
    throw FatalError( "main" ) << "No config file provided!";

  Debugging( "main" ) << "Reading config file stored in \"" << argv[1] << "\"";
  CepGen::Cards::PythonHandler card( argv[1] );
  mg.setParameters( card.parameters() );

  // We might want to cross-check visually the validity of our run
  mg.parameters->dump();

  // Let there be cross-section...
  double xsec, err;
  mg.computeXsection( xsec, err );

  CepGen::OutputHandler::HepMCHandler writer( "example.dat" );
  writer.setCrossSection( xsec, err );
  writer.initialise( *mg.parameters );

  // The events generation starts here !
  for ( unsigned int i = 0; i < mg.parameters->generation.maxgen; ++i ) {
    if ( i%1000 == 0 )
      cout << "Generating event #" << i+1 << endl;
    try {
      writer << mg.generateOneEvent().get();
    } catch ( CepGen::Exception& e ) { e.dump(); }
  }

  return 0;
}

