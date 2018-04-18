#ifndef CepGen_Core_Exception_h
#define CepGen_Core_Exception_h

#include <sstream>
#include <stdexcept>

#include "Logger.h"

#define CEPGEN_EXCEPT_MATCH( str ) \
  CepGen::Logger::get().passExceptionRule( str )

#define PrintMessage( mod ) \
  ( CepGen::Logger::get().level < CepGen::Logger::Nothing ) \
  ? CepGen::NullStream( mod ) \
  : CepGen::Exception( __PRETTY_FUNCTION__, mod, CepGen::kVerbatim )
#define Information( mod ) \
  ( CepGen::Logger::get().level < CepGen::Logger::Information && !CEPGEN_EXCEPT_MATCH( mod ) ) \
  ? CepGen::NullStream( mod ) \
  : CepGen::Exception( __PRETTY_FUNCTION__, mod, CepGen::kInformation )
#define Debugging( mod ) \
  ( CepGen::Logger::get().level < CepGen::Logger::Debug && !CEPGEN_EXCEPT_MATCH( mod ) ) \
  ? CepGen::NullStream( mod ) \
  : CepGen::Exception( __PRETTY_FUNCTION__, mod, CepGen::kDebugMessage )
#define DebuggingInsideLoop( mod ) \
  ( CepGen::Logger::get().level < CepGen::Logger::DebugInsideLoop && !CEPGEN_EXCEPT_MATCH( mod ) ) \
  ? CepGen::NullStream( mod ) \
  : CepGen::Exception( __PRETTY_FUNCTION__, mod, CepGen::kDebugMessage )
#define InWarning( mod ) \
  ( CepGen::Logger::get().level < CepGen::Logger::Warning && !CEPGEN_EXCEPT_MATCH( mod ) ) \
  ? CepGen::NullStream( mod ) \
  : CepGen::Exception( __PRETTY_FUNCTION__, mod, CepGen::kJustWarning )
#define InError( mod ) \
  ( CepGen::Logger::get().level < CepGen::Logger::Error && !CEPGEN_EXCEPT_MATCH( mod ) ) \
  ? CepGen::NullStream( mod ) \
  : CepGen::Exception( __PRETTY_FUNCTION__, mod, CepGen::kErrorMessage )
#define FatalError( mod ) \
  CepGen::Exception( __PRETTY_FUNCTION__, mod, CepGen::kErrorMessage )

namespace CepGen
{
  /// Enumeration of exception severities
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date 27 Mar 2015
  enum ExceptionType { kUndefined = -1, kDebugMessage, kVerbatim, kInformation, kJustWarning, kErrorMessage, kFatalError };

  /// A simple exception handler
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date 24 Mar 2015
  class Exception : public std::exception
  {
    public:
      /// Generic constructor
      /// \param[in] from method invoking the exception
      /// \param[in] module exception classifier
      /// \param[in] type exception type
      /// \param[in] id exception code (useful for logging)
      explicit inline Exception( const char* module = "", ExceptionType type = kUndefined, const int id = 0 ) :
        std::exception(), module_( module ), type_( type ), error_num_( id ) {}
      /// Generic constructor
      /// \param[in] from method invoking the exception
      /// \param[in] module exception classifier
      /// \param[in] type exception type
      /// \param[in] id exception code (useful for logging)
      explicit inline Exception( const char* from, const char* module, ExceptionType type = kUndefined, const int id = 0 ) :
        std::exception(), from_( from ), module_( module ), type_( type ), error_num_( id ) {}
      /// Generic constructor
      /// \param[in] from method invoking the exception
      /// \param[in] module exception classifier
      /// \param[in] type exception type
      /// \param[in] id exception code (useful for logging)
      explicit inline Exception( const char* from, const std::string& module, ExceptionType type = kUndefined, const int id = 0 ) :
        std::exception(), from_( from ), module_( module ), type_( type ), error_num_( id ) {}
      /// Copy constructor
      inline Exception( const Exception& rhs ) :
        from_( rhs.from_ ), module_( rhs.module_ ), type_( rhs.type_ ), error_num_( rhs.error_num_ ) {
        message_ << rhs.message_.str();
      }
      /// Default destructor (potentially killing the process)
      inline ~Exception() noexcept override {
        dump();
        // we stop this process' execution on fatal exception
        if ( type_ == kFatalError )
          exit(0);
      }

      //----- Overloaded stream operators

      /// Generic templated pipe operator
      template<typename T>
      inline friend const Exception& operator<<( const Exception& exc, T var ) {
        Exception& nc_except = const_cast<Exception&>( exc );
        nc_except.message_ << var;
        return exc;
      }
      /// Pipe modifier operator
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
          case kJustWarning: return "\033[34;1mJustWarning\033[0m";
          case kInformation: return "\033[32;1mInfo.\033[0m";
          case kDebugMessage: return "\033[33;1mDebug\033[0m";
          case kErrorMessage: return "\033[31;1mError\033[0m";
          case kFatalError: return "\033[31;1mFatal\033[0m";
          case kUndefined: default: return "\33[7;1mUndefined\033[0m";
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
        if ( type_ == kDebugMessage )
          os << " \033[30;4m" << from_ << "\033[0m\n";
        os << "\t" << message_.str();
        return os.str();
      }

    private:
      /// Extract a full exception message
      inline std::string fullMessage() const {
        if ( type_ == kInformation || type_ == kDebugMessage )
          return shortMessage();
        if ( type_ == kVerbatim )
          return message_.str();
        std::ostringstream os;
        os << "============================= Exception detected! =============================" << std::endl
           << " Class:       " << typeString() << std::endl;
        if ( !from_.empty() )
          os << " Raised by:   " << from_ << std::endl;
        os << " Description: \t" << message_.str() << std::endl;
        if ( errorNumber() != 0 )
          os << "-------------------------------------------------------------------------------" << std::endl
             << " Error #" << error_num_ << std::endl;
        os << "===============================================================================";
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
  /// Placeholder for debugging messages if logging threshold is not reached
  /// \date Apr 2018
  struct NullStream
  {
    explicit NullStream( const char* ) {}
    explicit NullStream( const std::string& ) {}
    NullStream( const Exception& ) {}
    template<class T> NullStream& operator<<( const T& ) { return *this; }
  };
}

#endif

