#include "CepGen/Processes/ProcessesHandler.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Core/ParametersList.h"
#include <iostream>

using namespace std;

int main( int argc, const char* argv[] )
{
  const char* proc_name = ( argc > 1 ) ? argv[1] : "lpair";
  cout << "Will build a process named \"" << proc_name << "\"." << endl;

  cout << "List of modules registered in the database:";
  for ( const auto& mod : cepgen::proc::ProcessesHandler::get().modules() )
    cout << " " << mod;
  cout << endl;

  auto proc = cepgen::proc::ProcessesHandler::get().build( proc_name, cepgen::ParametersList() );
  //--- at this point, the process has been found
  std::cout << "Successfully built the process \"" << proc->name() << "\"!\n"
    << " *) description: " << proc->description() << "\n"
    << " *) has event? " << proc->hasEvent() << "\n";
  if ( proc->hasEvent() ) { //--- dump a typical event content
    std::cout << "    event content (invalid kinematics, only check the parentage):\n";
    proc->addEventContent();
    proc->event()->dump();
  }

  return 0;
}
