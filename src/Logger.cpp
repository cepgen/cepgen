#include "Logger.h"
  
bool Logger::fBuilt = false;
Logger* Logger::fLogger = 0;

Logger*
Logger::GetInstance()
{
  if (!fBuilt) {
    fLogger = new Logger;
    fLogger->Level = Warning;
    fBuilt = true;
    return fLogger;
  }
  return fLogger;
}

