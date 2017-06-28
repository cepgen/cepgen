#ifndef Exception_h
#define Exception_h

#include <iomanip>
#include <sstream>

#include "Logger.h"

#define Information( m ) \
  if ( Logger::GetInstance()->Level>Logger::Nothing ) { Exception( __PRETTY_FUNCTION__, m, Information ).dump( Logger::GetInstance()->OutputStream ); }
#define Debugging( m ) \
  if ( Logger::GetInstance()->Level>=Logger::Debug )  { Exception( __PRETTY_FUNCTION__, m, DebugMessage ).dump( Logger::GetInstance()->OutputStream ); }
#define DebuggingInsideLoop( m ) \
  if ( Logger::GetInstance()->Level>=Logger::DebugInsideLoop ) { Exception( __PRETTY_FUNCTION__, m, DebugMessage ).dump( Logger::GetInstance()->OutputStream ); }
#define InWarning( m ) \
  if ( Logger::GetInstance()->Level>=Logger::Warning )  { Exception( __PRETTY_FUNCTION__, m, JustWarning ).dump( Logger::GetInstance()->OutputStream ); }
#define InError( m ) \
  if ( Logger::GetInstance()->Level>=Logger::Error )  { Exception( __PRETTY_FUNCTION__, m, ErrorMessage ).dump( Logger::GetInstance()->OutputStream ); }

/**
 * \brief Enumeration of exception severities
 * \author Laurent Forthomme <laurent.forthomme@cern.ch>
 * \date 27 Mar 2015
 */
typedef enum
{
  Undefined=-1,
  Information,
  DebugMessage,
  JustWarning,
  ErrorMessage,
  FatalError
} ExceptionType;

/**
 * \brief A simple exception handler
 * \author Laurent Forthomme <laurent.forthomme@cern.ch>
 * \date 24 Mar 2015
 */
class Exception
{
  public:
    /// Initialize a new exception object
    /// \param[in] from method invoking the exception
    /// \param[in] desc brief description of the exception
    /// \param[in] type exception type
    /// \param[in] id exception code (useful for logging)
    inline Exception( const char* from, std::string desc, ExceptionType type=Undefined, const int id=0 ) :
      from_( from ), description_( desc ), type_( type ), error_num_( id ) {}
    
    /// Initialize a new exception object
    /// \param[in] from method invoking the exception
    /// \param[in] desc brief description of the exception
    /// \param[in] type exception type
    /// \param[in] id exception code (useful for logging)
    inline Exception( const char* from, const char* desc, ExceptionType type=Undefined, const int id=0 ) :
      from_( from ), description_( desc ), type_( type ), error_num_( id ) {}

    inline ~Exception() {
      if ( type()==FatalError ) exit(0);
      // we stop this process' execution on fatal exception
    }
   
    /// Extract the origin of the exception 
    inline std::string from() const { return from_; }
    /// Extract the exception code
    inline int errorNumber() const { return error_num_; }
    /// Extract the brief exception description
    inline std::string description() const { return description_; }
    /// Extract the exception type
    inline ExceptionType type() const { return type_; }
    /// Extract a human-readable (and colourified) version of the exception type
    inline std::string typeString() const {
      switch ( type() ) {
        case JustWarning: return "\033[34;1mJustWarning\033[0m";
        case Information: return "\033[32;1mInfo\033[0m";
        case DebugMessage: return "\033[33;1mDebug\033[0m";
        case ErrorMessage: return "\033[31;1mError\033[0m";
        case FatalError: return "\033[31;1mFatal\033[0m";
        case Undefined: default: return "\33[7;1mUndefined\033[0m";
      }
    }
    
    /// Dump the full exception information in a given output stream
    /// \param[inout] os the output stream where the information is dumped
    inline void dump(std::ostream& os=std::cerr) const {
      if ( type()==Information ) {
        os << "[\033[32;1mInformation\033[0m]";
      }
      else if ( type()==DebugMessage ) {
        os << "==================================== \033[33;1mDebug\033[0m ====================================" << std::endl
           << " From:        " << from() << std::endl;
      }
      else {
        os << "============================= Exception detected! =============================" << std::endl
           << " Class:       " << typeString() << std::endl
           << " Raised by:   " << from() << std::endl;
      }
      if ( type()!=Information ) os << " Description: " << std::endl;
      os << "\t" << description() << std::endl;
      if ( errorNumber()!=0 )
        os << "-------------------------------------------------------------------------------" << std::endl
           << " Error #" << errorNumber() << std::endl;
      if ( type()!=Information ) os << "===============================================================================" << std::endl;
    }
    /// Extract a one-line summary of the exception
    inline std::string OneLine() const {
      std::ostringstream os;
      os << "[" << type() << "] === " << from() << " === "
         << description();
      return os.str();
    }
    
  private:
    /// Origin of the exception
    std::string from_;
    /// Small human-readable description
    std::string description_;
    /// Exception type
    ExceptionType type_;
    /// Integer exception number
    int error_num_;
};

#endif

