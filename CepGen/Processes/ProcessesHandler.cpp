#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"
#include "CepGen/Processes/PPtoFF.h"
#include "CepGen/Processes/ProcessesHandler.h"
#include <sstream>

namespace CepGen
{
  ProcessesHandler&
  ProcessesHandler::get()
  {
    static ProcessesHandler instance;
    return instance;
  }

  void
  ProcessesHandler::registerProcess( const std::string& name, const Process::GenericProcess* proc )
  {
    map_[name].reset( proc );
    CG_DEBUG( "ProcessesHandler" ) << "Process name \"" << name << "\" registered in database.";
  }

  void
  ProcessesHandler::dump() const
  {
    std::ostringstream oss;
    for ( const auto& p : map_ )
      oss << " '" << p.first << "'";
    CG_INFO( "ProcessesHandler:dump" )
      << "List of process(es) handled in the database:" << oss.str();
  }

  std::unique_ptr<Process::GenericProcess>
  ProcessesHandler::build( const std::string& name ) const
  {
    if ( map_.count( name ) == 0 )
      throw CG_FATAL( "ProcessesHandler:build" )
        << "Failed to retrieve a process with name \"" << name << "\"!";
    return map_.at( name )->clone();
  }
}
