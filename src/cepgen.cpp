#include "CepGen/Generator.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/CardsHandlerFactory.h"
#include "CepGen/Cards/Handler.h"

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

  cepgen::ArgumentsParser parser( argc, argv );
  parser
    .addArgument( "", "path to the configuration file", &input_card, 'i' )
    .addOptionalArgument( "num-events", "number of events to generate", -1, &num_events, 'n' )
    .addOptionalArgument( "list-modules", "list all runtime modules", false, &list_mods, 'l' )
    .addOptionalArgument( "debug", "debugging mode", false, &debug, 'd' )
    .parse();

  //--- first start by defining the generator object
  cepgen::Generator gen;

  //--- if modules listing is requested
  if ( list_mods ) {
    gen.dumpModules();
    return 0;
  }
  if ( debug )
    cepgen::utils::Logger::get().level = cepgen::utils::Logger::Level::debug;

  auto& params = gen.parameters();
  //--- no steering card nor additional flags found
  if ( input_card.empty() && parser.extra_config().empty() )
    throw CG_FATAL( "main" ) << "Neither input card nor configuration word provided!";
  //--- parse the steering card
  else if ( !input_card.empty() )
    params = cepgen::card::Handler::parse( input_card );
  //--- parse the additional flags
  else if ( !parser.extra_config().empty() )
    params = cepgen::card::CardsHandlerFactory::get().build( "cmd",
      cepgen::ParametersList().set<std::vector<std::string> >( "args", parser.extra_config() ) )
      ->parse( "", params );

  cepgen::utils::AbortHandler ctrl_c;

  try {
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
