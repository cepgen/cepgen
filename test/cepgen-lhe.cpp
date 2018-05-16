#include "CepGen/Cards/PythonHandler.h"
#include "CepGen/Cards/LpairHandler.h"

#include "CepGen/IO/HepMCHandler.h"
#include "CepGen/IO/LHEFHandler.h"

#include "CepGen/Generator.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"

#include <iostream>

using namespace std;

// we use polymorphism here
std::shared_ptr<CepGen::OutputHandler::ExportHandler> writer;

void storeEvent( const CepGen::Event& ev, unsigned long )
{
  if ( !writer )
    throw CG_FATAL( "storeEvent" ) << "Failed to retrieve a valid writer!";
  *writer << ev;
}

/**
 * Main caller for this Monte Carlo generator. Loads the configuration files'
 * variables if set as an argument to this program, else loads a default
 * "LHC-like" configuration, then launches the cross-section computation and
 * the events generation.
 * \author Laurent Forthomme <laurent.forthomme@cern.ch>
 */
int main( int argc, char* argv[] ) {

  if ( argc < 2 )
    throw CG_FATAL( "main" )
      << "No config file provided!\n\t"
      << "Usage: " << argv[0] << " config-file [format=lhef,hepmc] [filename=example.dat]";

  CepGen::Generator mg;

  //-----------------------------------------------------------------------------------------------
  // Steering card readout
  //-----------------------------------------------------------------------------------------------

  CG_DEBUG( "main" ) << "Reading config file stored in \"" << argv[1] << "\"";
  const string extension = CepGen::Cards::Handler::getExtension( argv[1] );
  if ( extension == "card" )
    mg.setParameters( CepGen::Cards::LpairHandler( argv[1] ).parameters() );
  else if ( extension == "py" )
    mg.setParameters( CepGen::Cards::PythonHandler( argv[1] ).parameters() );
  else
    throw CG_FATAL( "main" ) << "Unrecognized card format: ." << extension;

  //-----------------------------------------------------------------------------------------------
  // Output file writer definition
  //-----------------------------------------------------------------------------------------------

  const string format = ( argc > 2 ) ? argv[2] : "lhef";
  const char* filename = ( argc > 3 ) ? argv[3] : "example.dat";
  if ( format == "lhef" )
    writer = std::make_shared<CepGen::OutputHandler::LHEFHandler>( filename );
  else if ( format == "hepmc" )
    writer = std::make_shared<CepGen::OutputHandler::HepMCHandler>( filename );
  else
    throw CG_FATAL( "main" ) << "Unrecognized output format: " << format;

  //-----------------------------------------------------------------------------------------------
  // CepGen run part
  //-----------------------------------------------------------------------------------------------

  // We might want to cross-check visually the validity of our run
  mg.parameters->dump();

  // Let there be cross-section...
  double xsec = 0., err = 0.;
  mg.computeXsection( xsec, err );

  writer->initialise( *mg.parameters );
  writer->setCrossSection( xsec, err );

  // The events generation starts here!
  mg.generate( storeEvent );

  return 0;
}

