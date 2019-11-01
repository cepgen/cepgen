#include "CepGen/Cards/CardsHandler.h"
#include "CepGen/Processes/ProcessesHandler.h"
#include "CepGen/StructureFunctions/StructureFunctions.h"
#include "CepGen/Core/EventModifierHandler.h"
#include "CepGen/Core/ExportModuleHandler.h"

#include "CepGen/Generator.h"
#include "CepGen/Version.h"

#include <iostream>

using namespace std;
using namespace cepgen;

int main( int argc, const char* argv[] )
{
  cepgen::Generator gen;

  cout
    << "============================================================\n"
    << "CepGen v" << version() << "\n"
    << "List of modules registered in the runtime database:\n";
  {
    cout
      << "============================================================\n"
      << "Steering cards parsers definitions\n"
      << "------------------------------------------------------------\n";
    if ( card::CardsHandler::get().modules().empty() )
      cout << ">>> none found <<<" << endl;
    for ( const auto& mod : card::CardsHandler::get().modules() )
      cout << mod << "\n";
  }
  {
    cout
      << "------------------------------------------------------------\n"
      << "Processes definitions\n"
      << "------------------------------------------------------------\n";
    if ( proc::ProcessesHandler::get().modules().empty() )
      cout << ">>> none found <<<" << endl;
    for ( const auto& mod : proc::ProcessesHandler::get().modules() )
      cout << mod << " > " << proc::ProcessesHandler::get().build( mod )->description() << "\n";
  }
  {
    cout
      << "------------------------------------------------------------\n"
      << "Structure functions definitions\n"
      << "------------------------------------------------------------\n";
    if ( strfun::StructureFunctionsHandler::get().modules().empty() )
      cout << ">>> none found <<<" << endl;
    for ( const auto& mod : strfun::StructureFunctionsHandler::get().modules() )
      cout << mod << " > " << (strfun::Type)mod << "\n";
  }
  {
    cout
      << "------------------------------------------------------------\n"
      << "Event modification modules definitions\n"
      << "------------------------------------------------------------\n";
    if ( EventModifierHandler::get().modules().empty() )
      cout << ">>> none found <<<" << endl;
    for ( const auto& mod : EventModifierHandler::get().modules() )
      cout << mod << "\n";
  }
  {
    cout
      << "------------------------------------------------------------\n"
      << "Export modules definitions\n"
      << "------------------------------------------------------------\n";
    if ( io::ExportModuleHandler::get().modules().empty() )
      cout << ">>> none found <<<" << endl;
    for ( const auto& mod : io::ExportModuleHandler::get().modules() )
      cout << mod << "\n";
  }
  cout
    << "============================================================\n";

  return 0;
}
