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
  bool list_mods = false, cmd_line = false, debug = false;

  const std::string first_arg( argv[1] );
  if ( first_arg == "--cmd" || first_arg == "-c" )
    cmd_line = true;

  if ( !cmd_line )
    cepgen::ArgumentsParser( argc, argv )
      .addArgument( "", "path to the configuration file", &input_card, 'i' )
      .addOptionalArgument( "num-events", "number of events to generate", -1, &num_events, 'n' )
      .addOptionalArgument( "list-modules", "list all runtime modules", false, &list_mods, 'l' )
      .addOptionalArgument( "debug", "debugging mode", false, &debug, 'd' )
      .parse();

  if ( debug )
    cepgen::utils::Logger::get().level = cepgen::utils::Logger::Level::debug;

  cepgen::utils::AbortHandler ctrl_c;

  //--- first start by defining the generator object
  cepgen::Generator gen;

  if ( list_mods ) {
    gen.dumpModules();
    return 0;
  }

  try {
    //--- parse the steering card
    if ( cmd_line )
      gen.setParameters( cepgen::card::CardsHandlerFactory::get().build( "cmd",
        cepgen::ParametersList().set<std::vector<std::string> >( "args",
          std::vector<std::string>( argv+1, argv+argc ) ) )
        ->parameters() );
    else
      gen.setParameters( cepgen::card::Handler::parse( input_card )->parameters() );

    if ( num_events >= 0 ) { // user specified a number of events to generate
      gen.parameters().generation().maxgen = num_events;
      gen.parameters().generation().enabled = num_events > 0;
    }

    //--- list all parameters
    CG_LOG( "main" ) << gen.parametersPtr();

    //--- let there be a cross-section...
    double xsec = 0., err = 0.;
    gen.computeXsection( xsec, err );

    if ( gen.parameters().generation().enabled )
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
