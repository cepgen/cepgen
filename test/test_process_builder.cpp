#include "CepGen/Processes/ProcessesHandler.h"
#include "CepGen/Event/Event.h"
#include <iostream>

using namespace std;

int main( int argc, const char* argv[] )
{
  const char* proc_name = ( argc > 1 ) ? argv[1] : "lpair";
  cout << "Will build a process named \"" << proc_name << "\"." << endl;

  CepGen::ProcessesHandler::get().dump();

  auto proc = CepGen::ProcessesHandler::get().build( proc_name );
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
