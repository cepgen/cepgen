#ifndef CepGen_Core_Exception_h
#define CepGen_Core_Exception_h

#include <iomanip>
#include <sstream>

#include "Logger.h"

#define PrintMessage( m ) \
  if ( CepGen::Logger::get().level > CepGen::Logger::Nothing ) { CepGen::Exception( __PRETTY_FUNCTION__, m, CepGen::Verbatim ).dump( CepGen::Logger::get().outputStream ); }
#define Information( m ) \
  if ( CepGen::Logger::get().level > CepGen::Logger::Nothing ) { CepGen::Exception( __PRETTY_FUNCTION__, m, CepGen::Information ).dump( CepGen::Logger::get().outputStream ); }
#define Debugging( m ) \
  if ( CepGen::Logger::get().level >= CepGen::Logger::Debug )  { CepGen::Exception( __PRETTY_FUNCTION__, m, CepGen::DebugMessage ).dump( CepGen::Logger::get().outputStream ); }
#define DebuggingInsideLoop( m ) \
  if ( CepGen::Logger::get().level >= CepGen::Logger::DebugInsideLoop ) { CepGen::Exception( __PRETTY_FUNCTION__, m, CepGen::DebugMessage ).dump( CepGen::Logger::get().outputStream ); }
#define InWarning( m ) \
  if ( CepGen::Logger::get().level >= CepGen::Logger::Warning )  { CepGen::Exception( __PRETTY_FUNCTION__, m, CepGen::JustWarning ).dump( CepGen::Logger::get().outputStream ); }
#define InError( m ) \
  if ( CepGen::Logger::get().level >= CepGen::Logger::Error )  { CepGen::Exception( __PRETTY_FUNCTION__, m, CepGen::ErrorMessage ).dump( CepGen::Logger::get().outputStream ); }
#define FatalError( m ) \
  CepGen::Exception( __PRETTY_FUNCTION__, m, CepGen::FatalError ).dump( CepGen::Logger::get().outputStream );

namespace CepGen
{
  /**
   * \brief Enumeration of exception severities
   * \author Laurent Forthomme <laurent.forthomme@cern.ch>
   * \date 27 Mar 2015
   */
  enum ExceptionType { Undefined=-1, Verbatim, Information, DebugMessage, JustWarning, ErrorMessage, FatalError };

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
        if ( type() == FatalError ) exit(0);
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
        if ( type() == Verbatim ) {
          os << description() << std::endl;
          return;
        }
        if ( type() == Information ) {
          os << "[\033[32;1mInfo.\033[0m]\t" << description() << std::endl;
          return;
        }
        if ( type() == DebugMessage ) {
          os << "==================================== \033[33;1mDebug\033[0m ====================================" << std::endl
             << " From:        " << from() << std::endl;
        }
        else {
          os << "============================= Exception detected! =============================" << std::endl
             << " Class:       " << typeString() << std::endl
             << " Raised by:   " << from() << std::endl;
        }
        os << " Description: \t" << description() << std::endl;
        if ( errorNumber() != 0 )
          os << "-------------------------------------------------------------------------------" << std::endl
             << " Error #" << errorNumber() << std::endl;
        os << "===============================================================================" << std::endl;
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

  class Printer
  {
    public:
      inline static Exception LogInfo( const char* name ) { return Exception( "", name, Information ); }
      void operator<<( const char* text ) {
        LogInfo( text ).dump( Logger::get().outputStream );
      }
  };
}

#endif

