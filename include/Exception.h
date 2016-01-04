#ifndef Exception_h
#define Exception_h

#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib> // exit()

#include "Logger.h"

#define Information(m) \
  if (Logger::GetInstance()->Level>Logger::Nothing) { Exception(__PRETTY_FUNCTION__, m, Information).Dump(Logger::GetInstance()->OutputStream); }
#define Debug(m) \
  if (Logger::GetInstance()->Level>=Logger::Debug)  { Exception(__PRETTY_FUNCTION__, m, Debugging).Dump(Logger::GetInstance()->OutputStream); }
#define DebugInsideLoop(m) \
  if (Logger::GetInstance()->Level>=Logger::DebugInsideLoop) { Exception(__PRETTY_FUNCTION__, m, Debugging).Dump(Logger::GetInstance()->OutputStream); }
#define Warning(m) \
  if (Logger::GetInstance()->Level>=Logger::Warning)  { Exception(__PRETTY_FUNCTION__, m, JustWarning).Dump(Logger::GetInstance()->OutputStream); }
#define Error(m) \
  if (Logger::GetInstance()->Level>=Logger::Error)  { Exception(__PRETTY_FUNCTION__, m, Error).Dump(Logger::GetInstance()->OutputStream); }

/**
 * \brief Enumeration of exception severities
 * \author Laurent Forthomme <laurent.forthomme@cern.ch>
 * \date 27 Mar 2015
 */
typedef enum
{
  Undefined=-1,
  Information,
  Debugging,
  JustWarning,
  Error,
  Fatal
} ExceptionType;

/**
 * \brief A simple exception handler
 * \author Laurent Forthomme <laurent.forthomme@cern.ch>
 * \date 24 Mar 2015
 */
class Exception
{
  public:
    inline Exception(const char* from, std::string desc, ExceptionType type=Undefined, const int id=0) {
      fFrom = from;
      fDescription = desc;
      fType = type;
      fErrorNumber = id;
    }
    
    inline Exception(const char* from, const char* desc, ExceptionType type=Undefined, const int id=0) {
      fFrom = from;
      fDescription = desc;
      fType = type;
      fErrorNumber = id;
    }

    inline ~Exception() {
      if (Type()==Fatal) exit(0);
      // we stop this process' execution on fatal exception
    }
    
    inline std::string From() const { return fFrom; }
    inline int ErrorNumber() const { return fErrorNumber; }
    inline std::string Description() const { return fDescription; }
    inline ExceptionType Type() const { return fType; }
    inline std::string TypeString() const {
      switch (Type()) {
        case JustWarning: return "\033[34;1mJustWarning\033[0m";
        case Information: return "\033[33;1mInfo\033[0m";
        case Debugging: return "\033[32;1mDebug\033[0m";
        case Error: return "\033[31;1mError\033[0m";
        case Fatal: return "\033[31;1mFatal\033[0m";
        case Undefined: default: return "\33[7;1mUndefined\033[0m";
      }
    }
    
    inline void Dump(std::ostream& os=std::cerr) const {
      if (Type()==Information) {
        os << "================================= \033[33;1mInformation\033[0m =================================" << std::endl
           << " From:        " << From() << std::endl;
      }
      else if (Type()==Debugging) {
        os << "==================================== \033[32;1mDebug\033[0m ====================================" << std::endl
           << " From:        " << From() << std::endl;
      }
      else {
        os << "============================= Exception detected! =============================" << std::endl
           << " Class:       " << TypeString() << std::endl
           << " Raised by:   " << From() << std::endl;
      }
      os << " Description: " << std::endl
         << "\t" << Description() << std::endl;
      if (ErrorNumber()!=0)
        os << "-------------------------------------------------------------------------------" << std::endl
           << " Error #" << ErrorNumber() << std::endl;
      os << "===============================================================================" << std::endl;
    }
    inline std::string OneLine() const {
      std::ostringstream os;
      os << "[" << Type() << "] === " << From() << " === "
         << Description();
      return os.str();
    }
    
  private:
    std::string fFrom;
    std::string fDescription;
    ExceptionType fType;
    int fErrorNumber;
};

#endif

