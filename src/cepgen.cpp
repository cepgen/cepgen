#include "CepGen/Generator.h"

#include "CepGen/Modules/CardsHandlerFactory.h"
#include "CepGen/Cards/Handler.h"

#include "CepGen/Core/Exception.h"

#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/AbortHandler.h"

using namespace std;

/** Example executable for CepGen
 * - loads the steering card variables into the environment,
 * - launches the cross-section computation and the events generation (if requested).
 * \author Laurent Forthomme <laurent.forthomme@cern.ch>
 */
int main( int argc, char* argv[] )
{
  std::string input_card;
  int num_events;
  bool list_mods = false, debug = false;
  vector<string> addons;

  cepgen::ArgumentsParser parser( argc, argv );
  parser
    .addOptionalArgument( "config,i", "path to the configuration file", &input_card )
    .addOptionalArgument( "num-events,n", "number of events to generate", &num_events, -1 )
    .addOptionalArgument( "list-modules,l", "list all runtime modules", &list_mods, false )
    .addOptionalArgument( "add-ons,a", "external runtime plugin", &addons )
    .addOptionalArgument( "debug,d", "debugging mode", &debug, false )
    .parse();

  //--- first start by defining the generator object
  for ( const auto& lib : addons )
    try { cepgen::loadLibrary( lib ); } catch ( const cepgen::Exception& e ) {
      e.dump();
    }

  cepgen::Generator gen;

  //--- if modules listing is requested
  if ( list_mods ) {
    cepgen::dumpModules();
    return 0;
  }
  if ( debug )
    cepgen::utils::Logger::get().level = cepgen::utils::Logger::Level::debug;

  //--- no steering card nor additional flags found
  if ( input_card.empty() && parser.extra_config().empty() )
    throw CG_FATAL( "main" )
      << "Neither input card nor configuration word provided!\n\n "
      << parser.help_message();
  else {
    //--- parse the steering card
    if ( !input_card.empty() )
      gen.setParameters( cepgen::card::Handler::parse( input_card ) );
    //--- parse the additional flags
    if ( !parser.extra_config().empty() )
      gen.setParameters( cepgen::card::CardsHandlerFactory::get().build( "cmd",
        cepgen::ParametersList().set<std::vector<std::string> >( "args", parser.extra_config() ) )
        ->parse( "", &gen.parameters() ) );
  }

  cepgen::utils::AbortHandler ctrl_c;

  try {
    auto& params = gen.parameters();
    if ( num_events >= 0 ) { // user specified a number of events to generate
      params.generation().maxgen = num_events;
      params.generation().enabled = num_events > 0;
    }

    //--- list all parameters
    CG_LOG( "main" ) << gen.parametersPtr();

    //--- let there be a cross-section...
    double xsec = 0., err = 0.;
    gen.computeXsection( xsec, err );

    if ( params.generation().enabled )
      //--- events generation starts here
      // (one may use a callback function)
      gen.generate();
  } catch ( const cepgen::utils::RunAbortedException& e ) {
    CG_DEBUG( "main" ) << "Run aborted!";
  } catch ( const cepgen::Exception& e ) {
    e.dump();
  } catch ( const std::exception& e ) {
    CG_FATAL( "main" ) << "Other exception caught!\n\t" << e.what();
  }

  return 0;
}
