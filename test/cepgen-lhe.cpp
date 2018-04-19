#include "CepGen/Cards/PythonHandler.h"
#include "CepGen/IO/LHEFHandler.h"

#include "CepGen/Generator.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"

#include <iostream>

using namespace std;

std::shared_ptr<CepGen::OutputHandler::ExportHandler> writer;

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
    throw CG_FATAL( "main" ) << "No config file provided!";

  CG_DEBUG( "main" ) << "Reading config file stored in \"" << argv[1] << "\"";
  CepGen::Cards::PythonHandler card( argv[1] );
  mg.setParameters( card.parameters() );

  // We might want to cross-check visually the validity of our run
  mg.parameters->dump();

  // Let there be cross-section...
  double xsec = 0., err = 0.;
  mg.computeXsection( xsec, err );

  writer = std::make_shared<CepGen::OutputHandler::LHEFHandler>( "example.dat" );
  writer->initialise( *mg.parameters );
  writer->setCrossSection( xsec, err );

  // The events generation starts here!
  //FIXME move to a callback function for a more efficient usage of MT capabilities!
  for ( unsigned int i = 0; i < mg.parameters->generation.maxgen; ++i ) {
    if ( i % 1000 == 0 )
      cout << "Generating event #" << i+1 << endl;
    try {
      auto evt = *mg.generateOneEvent();
      *writer << evt;
    } catch ( CepGen::Exception& e ) { e.dump(); }
  }

  return 0;
}

