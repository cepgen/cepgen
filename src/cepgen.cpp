#include "CepGen/Generator.h"

#include "CepGen/Cards/CardsHandler.h"
#include "CepGen/Processes/ProcessesHandler.h"
#include "CepGen/StructureFunctions/StructureFunctions.h"

#include "CepGen/Modules/EventModifierHandler.h"
#include "CepGen/Modules/ExportModuleHandler.h"

#include "CepGen/Core/Exception.h"

#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/AbortHandler.h"

using namespace std;

void list_modules();

/** Example executable for CepGen
 * - loads the steering card variables into the environment,
 * - launches the cross-section computation and the events generation (if requested).
 * \author Laurent Forthomme <laurent.forthomme@cern.ch>
 */
int main( int argc, char* argv[] )
{
  std::string input_card;
  int num_events;
  bool list_mods;

  cepgen::ArgumentsParser( argc, argv )
    .addArgument( "", "path to the configuration file", &input_card, 'i' )
    .addOptionalArgument( "num-events", "number of events to generate", -1, &num_events, 'n' )
    .addOptionalArgument( "list-modules", "list all runtime modules", false, &list_mods, 'l' )
    .parse();

  //--- first start by defining the generator object
  cepgen::Generator gen;

  if ( list_mods ) {
    list_modules();
    return 0;
  }

  gen.setParameters( cepgen::card::Handler::parse( input_card )->parameters() );

  if ( num_events >= 0 ) { // user specified a number of events to generate
    gen.parameters().generation().maxgen = num_events;
    gen.parameters().generation().enabled = num_events > 0;
  }

  //--- list all parameters
  CG_LOG( "main" ) << gen.parametersPtr();

  cepgen::utils::AbortHandler ctrl_c;

  try {
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
    CG_FATAL( "main" ) << "Other exception caught!\n\t"
      << e.what();
  }

  return 0;
}

void list_modules()
{
  using namespace cepgen;

  string sep_mid( 80, '-' ), sep_big( 80, '=' );
  sep_mid += "\n", sep_big += "\n";

  cout << sep_big
    << "List of modules registered in the runtime database:\n";
  {
    cout << sep_big << "Steering cards parsers definitions\n" << sep_mid;
    if ( cepgen::card::CardsHandler::get().modules().empty() )
      cout << ">>> none found <<<" << endl;
    for ( const auto& mod : cepgen::card::CardsHandler::get().modules() )
      cout << "." << mod << " extension\n";
  }
  {
    cout << sep_mid << "Processes definitions\n" << sep_mid;
    if ( cepgen::proc::ProcessesHandler::get().modules().empty() )
      cout << ">>> none found <<<" << endl;
    for ( const auto& mod : cepgen::proc::ProcessesHandler::get().modules() )
      cout << mod << " > " << cepgen::proc::ProcessesHandler::get().build( mod )->description() << "\n";
  }
  {
    cout << sep_mid << "Structure functions definitions\n" << sep_mid;
    if ( cepgen::strfun::StructureFunctionsHandler::get().modules().empty() )
      cout << ">>> none found <<<" << endl;
    for ( const auto& mod : cepgen::strfun::StructureFunctionsHandler::get().modules() )
      cout << mod << " > " << (cepgen::strfun::Type)mod << "\n";
  }
  {
    cout << sep_mid << "Event modification modules definitions\n" << sep_mid;
    if ( cepgen::EventModifierHandler::get().modules().empty() )
      cout << ">>> none found <<<" << endl;
    for ( const auto& mod : cepgen::EventModifierHandler::get().modules() )
      cout << mod << "\n";
  }
  {
    cout << sep_mid << "Export modules definitions\n" << sep_mid;
    if ( cepgen::io::ExportModuleHandler::get().modules().empty() )
      cout << ">>> none found <<<" << endl;
    for ( const auto& mod : cepgen::io::ExportModuleHandler::get().modules() )
      cout << mod << "\n";
  }
  cout << sep_big;
}
