#include "CepGen/Parameters.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/ArgumentsParser.h"

#include "CepGen/Cards/Handler.h"
#include "CepGen/Modules/CardsHandlerFactory.h"

using namespace std;

int main( int argc, char* argv[] )
{
  string input_config, output_config;

  cepgen::ArgumentsParser parser( argc, argv );
  parser
    .addArgument( "input", "input configuration", &input_config, 'i' )
    .addArgument( "output", "output output", &output_config, 'o' )
    .parse();

  try {
    auto params = cepgen::card::Handler::parse( input_config );
    cepgen::card::Handler::write( params, output_config );
    CG_INFO( "main" )
      << "Successfully converted the \""
      << cepgen::card::Handler::extension( input_config )
      << "\" card into a \""
      << cepgen::card::Handler::extension( output_config )
      << "\" card.\n\t"
      << "\"" << output_config << "\" file created.";

  } catch ( const cepgen::Exception& e ) {
    throw CG_FATAL( "main" )
      << "Failed to convert a \""
      << cepgen::card::Handler::extension( input_config )
      << "\" card into a \""
      << cepgen::card::Handler::extension( output_config )
      << "\" card!\n"
      << e.message();
  }

  return 0;
}
