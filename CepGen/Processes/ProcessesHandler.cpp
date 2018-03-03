#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"
#include "CepGen/Processes/PPtoLL.h"
#include "CepGen/Processes/ProcessesHandler.h"

namespace CepGen
{
  ProcessesHandler&
  ProcessesHandler::get()
  {
    static ProcessesHandler instance;
    return instance;
  }

  ProcessesHandler::ProcessesHandler()
  {
  }

  ProcessesHandler::~ProcessesHandler()
  {}

  CepGen::Process::GenericProcess*
  ProcessesHandler::build( const char* name ) const
  {
    if ( map_.count( name ) == 0 )
      return nullptr;
    return new Process::GenericProcess( map_.at( name ) );
  }
}
