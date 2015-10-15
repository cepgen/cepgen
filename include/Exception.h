#ifndef Exception_h
#define Exception_h

#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib> // exit()

enum LoggingLevel { Nothing, Error, Warning, Information, Debug, DebugInsideLoop };
//extern LoggingLevel kLoggingLevel;

#define Info(m)  if (kLoggingLevel>Nothing) { Exception(__PRETTY_FUNCTION__, m, Info).Dump(); }
#define Debug(m) if (kLoggingLevel>=Debug)  { Exception(__PRETTY_FUNCTION__, m, Debugging).Dump(); }
#define DebugInsideLoop(m) if (kLoggingLevel>=DebugInsideLoop) { Exception(__PRETTY_FUNCTION__, m, Debugging).Dump(); }

/**
 * \brief Enumeration of exception severities
 * \author Laurent Forthomme <laurent.forthomme@cern.ch>
 * \date 27 Mar 2015
 */
typedef enum
{
  Undefined=-1,
  Info,
  Debugging,
  JustWarning,
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
        case Info: return "\033[33;1mInfo\033[0m";
        case Debugging: return "\033[32;1mDebug\033[0m";
        case Fatal: return "\033[31;1mFatal\033[0m";
        case Undefined: default: return "\33[7;1mUndefined\033[0m";
      }
    }
    
    inline void Dump(std::ostream& os=std::cerr) const {
      if (Type()==Info) {
        os << "============================ \033[33;1mInformation\033[0m ============================" << std::endl
           << " From:        " << From() << std::endl;
      }
      else {
        os << "======================== Exception detected! ========================" << std::endl
           << " Class:       " << TypeString() << std::endl
           << " Raised by:   " << From() << std::endl;
      }
      os << " Description: " << std::endl
         << "\t" << Description() << std::endl;
      if (ErrorNumber()!=0)
        os << "---------------------------------------------------------------------" << std::endl
           << " Error #" << ErrorNumber() << std::endl;
      os << "=====================================================================" << std::endl;
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

static LoggingLevel kLoggingLevel = Warning;
