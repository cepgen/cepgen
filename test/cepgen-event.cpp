#include "CepGen/Cards/Handler.h"

#include "CepGen/IO/ExportHandler.h"

#include "CepGen/Generator.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Event/Event.h"

#include "AbortHandler.h"

#include <iostream>

using namespace std;

// we use polymorphism here
std::unique_ptr<cepgen::output::GenericExportHandler> writer;

void storeEvent( const cepgen::Event& ev, unsigned long )
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
int main( int argc, char* argv[] )
{
  ostringstream os;
  string delim;
  for ( const auto& mod : cepgen::output::ExportHandler::get().modules() )
    os << delim << mod, delim = ",";

  if ( argc < 2 )
    throw CG_FATAL( "main" )
      << "No config file provided!\n\t"
      << "Usage: " << argv[0] << " config-file [format=" << os.str() << "] [filename=example.dat]";

  cepgen::Generator mg;

  //-----------------------------------------------------------------------------------------------
  // Steering card readout
  //-----------------------------------------------------------------------------------------------

  CG_DEBUG( "main" ) << "Reading config file stored in \"" << argv[1] << "\"";
  mg.setParameters( cepgen::card::Handler::parse( argv[1] ) );

  //-----------------------------------------------------------------------------------------------
  // Output file writer definition
  //-----------------------------------------------------------------------------------------------

  const string format = ( argc > 2 ) ? argv[2] : "lhef";
  const char* filename = ( argc > 3 ) ? argv[3] : "example.dat";
  writer = cepgen::output::ExportHandler::get().build( format, cepgen::ParametersList()
    .set<std::string>( "filename", filename ) );

  //-----------------------------------------------------------------------------------------------
  // CepGen run part
  //-----------------------------------------------------------------------------------------------

  // We might want to cross-check visually the validity of our run
  CG_INFO( "main" ) << mg.parametersPtr();
  cepgen::utils::AbortHandler ctrl_c;

  //--- let there be cross-section...
  double xsec = 0., err = 0.;
  mg.computeXsection( xsec, err );

  writer->initialise( mg.parameters() );
  writer->setCrossSection( xsec, err );

  // The events generation starts here!
  mg.generate( storeEvent );

  return 0;
}
