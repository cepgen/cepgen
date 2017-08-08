#include "Logger.h"
  
bool Logger::built_ = false;
Logger* Logger::logger_ = 0;

Logger*
Logger::GetInstance()
{
  if ( built_ ) return logger_;

  logger_ = new Logger;
  logger_->Level = Error;
  built_ = true;

  return logger_;
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
