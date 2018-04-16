#ifndef CepGen_Core_Exception_h
#define CepGen_Core_Exception_h

#include <sstream>
#include <stdexcept>

#include "Logger.h"

/*#define PrintMessage( m ) \
  { if ( CepGen::Logger::get().level > CepGen::Logger::Nothing ) { CepGen::Exception( __PRETTY_FUNCTION__, m, CepGen::Verbatim ).dump( CepGen::Logger::get().outputStream ); } }
#define Information( m ) \
  { if ( CepGen::Logger::get().level > CepGen::Logger::Nothing ) { CepGen::Exception( __PRETTY_FUNCTION__, m, CepGen::Information ).dump( CepGen::Logger::get().outputStream ); } }
#define Debugging( m ) \
  { if ( CepGen::Logger::get().level >= CepGen::Logger::Debug )  { CepGen::Exception( __PRETTY_FUNCTION__, m, CepGen::DebugMessage ).dump( CepGen::Logger::get().outputStream ); } }
#define DebuggingInsideLoop( m ) \
  { if ( CepGen::Logger::get().level >= CepGen::Logger::DebugInsideLoop ) { CepGen::Exception( __PRETTY_FUNCTION__, m, CepGen::DebugMessage ).dump( CepGen::Logger::get().outputStream ); } }
#define InWarning( m ) \
  { if ( CepGen::Logger::get().level >= CepGen::Logger::Warning )  { CepGen::Exception( __PRETTY_FUNCTION__, m, CepGen::JustWarning ).dump( CepGen::Logger::get().outputStream ); } }
#define InError( m ) \
  { if ( CepGen::Logger::get().level >= CepGen::Logger::Error )  { CepGen::Exception( __PRETTY_FUNCTION__, m, CepGen::ErrorMessage ).dump( CepGen::Logger::get().outputStream ); } }
#define FatalError( m ) \
  { CepGen::Exception( __PRETTY_FUNCTION__, m, CepGen::FatalError ).dump( CepGen::Logger::get().outputStream ); }*/

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
  class Exception : public std::exception
  {
    public:
      /// Initialize a new exception object
      /// \param[in] from method invoking the exception
      /// \param[in] module exception classifier
      /// \param[in] type exception type
      /// \param[in] id exception code (useful for logging)
      explicit inline Exception( const char* module, ExceptionType type = Undefined, const int id = 0 ) :
        std::exception(), module_( module ), type_( type ), error_num_( id ) {}
      /// Initialize a new exception object
      /// \param[in] from method invoking the exception
      /// \param[in] module exception classifier
      /// \param[in] type exception type
      /// \param[in] id exception code (useful for logging)
      explicit inline Exception( const char* from, const char* module, ExceptionType type = Undefined, const int id = 0 ) :
        std::exception(), from_( from ), module_( module ), type_( type ), error_num_( id ) {}
      /// Initialize a new exception object
      /// \param[in] from method invoking the exception
      /// \param[in] module exception classifier
      /// \param[in] type exception type
      /// \param[in] id exception code (useful for logging)
      explicit inline Exception( const char* from, const std::string& module, ExceptionType type = Undefined, const int id = 0 ) :
        std::exception(), from_( from ), module_( module ), type_( type ), error_num_( id ) {}

      inline Exception( const Exception& rhs ) :
        from_( rhs.from_ ), module_( rhs.module_ ), type_( rhs.type_ ), error_num_( rhs.error_num_ ) {
        message_ << rhs.message_.str();
      }

      inline ~Exception() noexcept override {
        dump();
        if ( type() == FatalError )
          exit(0);
        // we stop this process' execution on fatal exception
      }

      //----- Overloaded stream operators

      template<typename T>
      inline friend const Exception& operator<<( const Exception& exc, T var ) {
        Exception& nc_except = const_cast<Exception&>( exc );
        nc_except.message_ << var;
        return exc;
      }
      inline friend const Exception& operator<<( const Exception& exc, std::ios_base&( *f )( std::ios_base& ) ) {
        Exception& nc_except = const_cast<Exception&>( exc );
        f( nc_except.message_ );
        return exc;
      }

      /// Exception message
      inline const char* what() const noexcept override {
        return fullMessage().c_str();
      }

      /// Extract the origin of the exception
      inline std::string from() const { return from_; }
      /// Extract the exception code
      inline int errorNumber() const { return error_num_; }
      /// Extract the exception type
      inline ExceptionType type() const { return type_; }
      /// Extract a human-readable (and colourified) version of the exception type
      inline std::string typeString() const {
        switch ( type() ) {
          case JustWarning: return "\033[34;1mJustWarning\033[0m";
          case Information: return "\033[32;1mInfo.\033[0m";
          case DebugMessage: return "\033[33;1mDebug\033[0m";
          case ErrorMessage: return "\033[31;1mError\033[0m";
          case FatalError: return "\033[31;1mFatal\033[0m";
          case Undefined: default: return "\33[7;1mUndefined\033[0m";
        }
      }

      /// Dump the full exception information in a given output stream
      /// \param[inout] os the output stream where the information is dumped
      inline void dump( std::ostream& os = Logger::get().outputStream ) {
        os << fullMessage() << std::endl;
      }
      /// Extract a one-line summary of the exception
      inline std::string shortMessage() const {
        std::ostringstream os;
        os << "[" << typeString() << "]";
        if ( type_ == DebugMessage )
          os << " \033[30;4m" << from_ << "\033[0m\n";
        os << "\t" << message_.str();
        return os.str();
      }

    private:
      inline std::string fullMessage() const {
        if ( type_ == Information || type_ == DebugMessage )
          return shortMessage();
        if ( type_ == Verbatim )
          return message_.str();
        std::ostringstream os;
        os << "============================= Exception detected! =============================" << std::endl
           << " Class:       " << typeString() << std::endl
           << " Raised by:   " << from_ << std::endl;
        os << " Description: \t" << message_.str() << std::endl;
        if ( errorNumber() != 0 )
          os << "-------------------------------------------------------------------------------" << std::endl
             << " Error #" << error_num_ << std::endl;
        os << "===============================================================================" << std::endl;
        return os.str();
      }
      /// Origin of the exception
      std::string from_;
      /// Exception classificator
      std::string module_;
      /// Message to throw
      std::ostringstream message_;
      /// Exception type
      ExceptionType type_;
      /// Integer exception number
      int error_num_;
  };

#define FatalError( mod ) CepGen::Exception( __PRETTY_FUNCTION__, mod, CepGen::FatalError )
#define InError( mod ) CepGen::Exception( __PRETTY_FUNCTION__, mod, CepGen::ErrorMessage )
#define InWarning( mod ) CepGen::Exception( __PRETTY_FUNCTION__, mod, CepGen::JustWarning )
#define Information( mod ) CepGen::Exception( __PRETTY_FUNCTION__, mod, CepGen::Information )
#define PrintMessage( mod ) CepGen::Exception( __PRETTY_FUNCTION__, mod, CepGen::Verbatim )
#define Debugging( mod ) CepGen::Exception( __PRETTY_FUNCTION__, mod, CepGen::DebugMessage )
#define DebuggingInsideLoop( mod ) CepGen::Exception( __PRETTY_FUNCTION__, mod, CepGen::DebugMessage )
}

#endif

