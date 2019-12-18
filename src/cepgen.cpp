#include "CepGen/Generator.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Cards/Handler.h"

#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/AbortHandler.h"

#include "CepGen/Modules/CardsHandlerFactory.h"

#include "CepGen/Modules/ProcessesFactory.h"
#include "CepGen/Modules/Process.h"

#include "CepGen/Physics/AlphaS.h"

#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/StructureFunctions/SigmaRatio.h"

#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Modules/ExportModuleFactory.h"

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

  cepgen::utils::AbortHandler ctrl_c;

  try {
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
  { cout << sep_big << "Steering cards parsers definitions\n" << sep_mid;
    if ( cepgen::card::CardsHandlerFactory::get().modules().empty() )
      cout << ">>> none found <<<" << endl;
    for ( const auto& mod : cepgen::card::CardsHandlerFactory::get().modules() )
      cout << "." << mod << " extension\n";
  }
  { cout << sep_mid << "Processes definitions\n" << sep_mid;
    if ( cepgen::proc::ProcessesFactory::get().modules().empty() )
      cout << ">>> none found <<<" << endl;
    for ( const auto& mod : cepgen::proc::ProcessesFactory::get().modules() )
      cout << mod << " > " << cepgen::proc::ProcessesFactory::get().build( mod )->description() << "\n";
  }
  { cout << sep_mid << "Structure functions definitions\n" << sep_mid;
    if ( cepgen::strfun::StructureFunctionsFactory::get().modules().empty() )
      cout << ">>> none found <<<" << endl;
    for ( const auto& mod : cepgen::strfun::StructureFunctionsFactory::get().modules() )
      cout << mod << " > " << (cepgen::strfun::Type)mod << "\n";
  }
  { cout << sep_mid << "Cross section ratios definitions\n" << sep_mid;
    if ( cepgen::sigrat::SigmaRatiosFactory::get().modules().empty() )
      cout << ">>> none found <<<" << endl;
    for ( const auto& mod : cepgen::sigrat::SigmaRatiosFactory::get().modules() )
      cout << mod << " > " << (cepgen::sigrat::Type)mod << "\n";
  }
  { cout << sep_mid << "Event modification modules definitions\n" << sep_mid;
    if ( cepgen::EventModifierFactory::get().modules().empty() )
      cout << ">>> none found <<<" << endl;
    for ( const auto& mod : cepgen::EventModifierFactory::get().modules() )
      cout << mod << "\n";
  }
  { cout << sep_mid << "Export modules definitions\n" << sep_mid;
    if ( cepgen::io::ExportModuleFactory::get().modules().empty() )
      cout << ">>> none found <<<" << endl;
    for ( const auto& mod : cepgen::io::ExportModuleFactory::get().modules() )
      cout << mod << "\n";
  }
  { cout << sep_mid << "alpha(s) evolution algorithms definitions\n" << sep_mid;
    if ( cepgen::AlphaSFactory::get().modules().empty() )
      cout << ">>> none found <<<" << endl;
    for ( const auto& mod : cepgen::AlphaSFactory::get().modules() )
      cout << mod << "\n";
  }
  cout << sep_big;
}
