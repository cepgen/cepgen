#include <iostream>

#include "CepGen/Generator.h"
#include "CepGen/Export/EventWriter.h"
#include "CepGen/Cards/LpairHandler.h"

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
  
  if ( argc==1 ) InError( "No config file provided." );

  Debugging( Form( "Reading config file stored in %s", argv[1] ) );
  CepGen::Cards::LpairHandler card( argv[1] );
  mg.setParameters( card.parameters() );

  // We might want to cross-check visually the validity of our run
  mg.parameters->dump();

  // Let there be cross-section...
  double xsec, err;
  mg.computeXsection( xsec, err );

  //if ( !mg.parameters->generation ) return 0;

  //CepGen::OutputHandler::EventWriter writer( CepGen::OutputHandler::ExportHandler::LHE, "example.dat" );
  CepGen::OutputHandler::EventWriter writer( CepGen::OutputHandler::ExportHandler::HepMC, "example.dat" );
  writer.setCrossSection( xsec, err );
  writer.initialise( *mg.parameters );

  Information( Form( "HepMC version: %s", HepMC::versionName().c_str() ) );

  // The events generation starts here !
  for ( unsigned int i=0; i<mg.parameters->generation.maxgen; i++ ) {
    if ( i%1000==0 )
      cout << "Generating event #" << i+1 << endl;
    try {
      const CepGen::Event* ev = mg.generateOneEvent();
      writer << ev;
    } catch ( CepGen::Exception& e ) { e.dump(); }
  }

  //mg.parameters->StoreConfigFile( "lastrun.card" );

  return 0;
}

