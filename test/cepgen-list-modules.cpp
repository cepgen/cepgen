#include "CepGen/Processes/ProcessesHandler.h"
#include "CepGen/StructureFunctions/StructureFunctions.h"
#include "CepGen/Core/EventModifierHandler.h"
#include "CepGen/IO/ExportHandler.h"

#include "CepGen/Version.h"

#include <iostream>

using namespace std;
using namespace cepgen;

int main( int argc, const char* argv[] )
{
  cout
    << "============================================================\n"
    << "CepGen v" << version() << "\n"
    << "List of modules registered in the runtime database:\n";
  {
    cout
      << "============================================================\n"
      << "Processes definitions\n"
      << "------------------------------------------------------------\n";
    for ( const auto& mod : proc::ProcessesHandler::get().modules() )
      cout << mod << " > " << proc::ProcessesHandler::get().build( mod )->description() << "\n";
  }
  {
    cout
      << "------------------------------------------------------------\n"
      << "Structure functions definitions\n"
      << "------------------------------------------------------------\n";
    for ( const auto& mod : strfun::StructureFunctionsHandler::get().modules() )
      cout << mod << " > " << (strfun::Type)mod << "\n";
  }
  {
    cout
      << "------------------------------------------------------------\n"
      << "Event modification modules definitions\n"
      << "------------------------------------------------------------\n";
    for ( const auto& mod : EventModifierHandler::get().modules() )
      cout << mod << "\n";
  }
  {
    cout
      << "------------------------------------------------------------\n"
      << "Export modules definitions\n"
      << "------------------------------------------------------------\n";
    for ( const auto& mod : io::ExportHandler::get().modules() )
      cout << mod << "\n";
  }
  cout
    << "============================================================\n";

  return 0;
}
