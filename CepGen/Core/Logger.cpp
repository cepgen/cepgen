#include "Logger.h"

namespace CepGen
{
  Logger&
  Logger::get()
  {
    static Logger log;
    return log;
  }

  std::ostream&
  operator<<( std::ostream& os, const Logger::LoggingLevel& lvl )
  {
    switch ( lvl ) {
      case Logger::Nothing:         os << "None"; break;
      case Logger::Error:           os << "Errors"; break;
      case Logger::Warning:         os << "Warnings"; break;
      case Logger::Information:     os << "Infos"; break;
      case Logger::Debug:           os << "Debug"; break;
      case Logger::DebugInsideLoop: os << "Debug (in loops)"; break;
    }
    return os;
  }
}
