#include "Logger.h"
  
bool Logger::fBuilt = false;
Logger* Logger::fLogger = 0;

Logger*
Logger::GetInstance()
{
  if (!fBuilt) {
    fLogger = new Logger;
    fLogger->Level = Error;
    fBuilt = true;
    return fLogger;
  }
  return fLogger;
}

std::ostream&
operator<<(std::ostream& os, const Logger::LoggingLevel& lvl)
{
  switch (lvl) {
  case Logger::Nothing:         os << "None"; break;
  case Logger::Error:           os << "Errors"; break;
  case Logger::Warning:         os << "Warnings"; break;
  case Logger::Information:     os << "Infos"; break;
  case Logger::Debug:           os << "Debug"; break;
  case Logger::DebugInsideLoop: os << "Debug (in loops)"; break;
  }
  return os;
}
