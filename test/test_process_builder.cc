#include "CepGen/Modules/ProcessesFactory.h"
#include "CepGen/Processes/Process.h"

#include "CepGen/Generator.h"
#include "CepGen/Core/ParametersList.h"

#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/String.h"

#include <iostream>

using namespace std;

int main( int argc, char* argv[] )
{
  string proc_name;
  bool list;

  cepgen::ArgumentsParser( argc, argv )
    .addOptionalArgument( "proc-name,p", "name of the process", &proc_name )
    .addOptionalArgument( "list,l", "list all processes", &list, false )
    .parse();

  cepgen::initialise();

  if ( list ) {
    cout << "List of modules registered in the runtime database:";
    for ( const auto& mod : cepgen::proc::ProcessesFactory::get().modules() )
      cout << "\n> " << cepgen::utils::boldify( mod );
    cout << endl;
    return 0;
  }

  if ( !proc_name.empty() ) {
    cout << "Will build a process named \"" << proc_name << "\"." << endl;

    auto proc = cepgen::proc::ProcessesFactory::get().build( proc_name, cepgen::ParametersList() );
    //--- at this point, the process has been found
    std::cout << "Successfully built the process \"" << proc->name() << "\"!\n"
      << " *) description: " << proc->description() << "\n"
      << " *) has event? " << proc->hasEvent() << "\n";
    if ( proc->hasEvent() ) { //--- dump a typical event content
      std::cout << "    event content (invalid kinematics, only check the parentage):\n";
      proc->addEventContent();
      proc->event().dump();
    }
  }

  return 0;
}
