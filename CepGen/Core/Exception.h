#ifndef CepGen_Core_Exception_h
#define CepGen_Core_Exception_h

#include <sstream>
#include <stdexcept>
#include <csignal>

#include "CepGen/Core/Logger.h"

#define CG_EXCEPT_MATCH( str, type ) \
  CepGen::Logger::get().passExceptionRule( str, CepGen::Logger::Level::type )

#define CG_LOG( mod ) \
  ( !CG_EXCEPT_MATCH( mod, information ) ) \
  ? CepGen::NullStream( mod ) \
  : CepGen::Exception( __PRETTY_FUNCTION__, mod, CepGen::Exception::Type::verbatim )
#define CG_INFO( mod ) \
  ( !CG_EXCEPT_MATCH( mod, information ) ) \
  ? CepGen::NullStream( mod ) \
  : CepGen::Exception( __PRETTY_FUNCTION__, mod, CepGen::Exception::Type::info )
#define CG_DEBUG( mod ) \
  ( !CG_EXCEPT_MATCH( mod, debug ) ) \
  ? CepGen::NullStream( mod ) \
  : CepGen::Exception( __PRETTY_FUNCTION__, mod, CepGen::Exception::Type::debug )
#define CG_DEBUG_LOOP( mod ) \
  ( !CG_EXCEPT_MATCH( mod, debugInsideLoop ) ) \
  ? CepGen::NullStream( mod ) \
  : CepGen::Exception( __PRETTY_FUNCTION__, mod, CepGen::Exception::Type::debug )
#define CG_WARNING( mod ) \
  ( !CG_EXCEPT_MATCH( mod, warning ) ) \
  ? CepGen::NullStream( mod ) \
  : CepGen::Exception( __PRETTY_FUNCTION__, mod, CepGen::Exception::Type::warning )
#define CG_ERROR( mod ) \
  ( !CG_EXCEPT_MATCH( mod, error ) ) \
  ? CepGen::NullStream( mod ) \
  : CepGen::Exception( __PRETTY_FUNCTION__, mod, CepGen::Exception::Type::warning )
#define CG_FATAL( mod ) \
  CepGen::Exception( __PRETTY_FUNCTION__, mod, CepGen::Exception::Type::fatal )

namespace CepGen
{
  /// \brief A simple exception handler
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date 24 Mar 2015
  class Exception : public std::exception
  {
    public:
      /// Enumeration of exception severities
      /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
      /// \date 27 Mar 2015
      enum class Type {
        undefined = -1, debug, verbatim, info, warning, error, fatal };

      /// Generic constructor
      /// \param[in] module exception classifier
      /// \param[in] type exception type
      /// \param[in] id exception code (useful for logging)
      explicit inline Exception( const char* module = "", Type type = Type::undefined, const int id = 0 ) :
        module_( module ), type_( type ), error_num_( id ) {}
      /// Generic constructor
      /// \param[in] from method invoking the exception
      /// \param[in] module exception classifier
      /// \param[in] type exception type
      /// \param[in] id exception code (useful for logging)
      explicit inline Exception( const char* from, const char* module, Type type = Type::undefined, const int id = 0 ) :
        from_( from ), module_( module ), type_( type ), error_num_( id ) {}
      /// Generic constructor
      /// \param[in] from method invoking the exception
      /// \param[in] module exception classifier
      /// \param[in] type exception type
      /// \param[in] id exception code (useful for logging)
      explicit inline Exception( const char* from, const std::string& module, Type type = Type::undefined, const int id = 0 ) :
        from_( from ), module_( module ), type_( type ), error_num_( id ) {}
      /// Copy constructor
      inline Exception( const Exception& rhs ) :
        from_( rhs.from_ ), module_( rhs.module_ ), message_( rhs.message_.str() ), type_( rhs.type_ ), error_num_( rhs.error_num_ ) {}
      /// Default destructor (potentially killing the process)
      inline ~Exception() noexcept override {
        do { dump(); } while ( 0 );
        // we stop this process' execution on fatal exception
        if ( type_ == Type::fatal )
          if ( raise( SIGINT ) != 0 )
            exit( 0 );
      }

      //----- Overloaded stream operators

      /// Generic templated message feeder operator
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
      /*inline const char* what() const noexcept override {
        return message_.str().c_str();
      }*/

      /// Extract the origin of the exception
      inline std::string from() const { return from_; }
      /// Extract the exception code
      inline int errorNumber() const { return error_num_; }
      /// Extract the exception type
      inline Type type() const { return type_; }
      /// Extract a human-readable (and colourified) version of the exception type
      inline std::string typeString() const {
        switch ( type() ) {
          case Type::warning: return "\033[34;1mWarning\033[0m";
          case Type::info: return "\033[32;1mInfo.\033[0m";
          case Type::debug: return "\033[33;1mDebug\033[0m";
          case Type::error: return "\033[31;1mError\033[0m";
          case Type::fatal: return "\033[31;1mFatal\033[0m";
          case Type::undefined: default: return "\33[7;1mUndefined\033[0m";
        }
      }

      /// Dump the full exception information in a given output stream
      /// \param[inout] os the output stream where the information is dumped
      inline void dump( std::ostream& os = *Logger::get().output ) const {
        os << fullMessage() << std::endl;
      }
      /// Extract a one-line summary of the exception
      inline std::string shortMessage() const {
        std::ostringstream os;
        os << "[" << typeString() << "]";
        if ( type_ == Type::debug || type_ == Type::warning )
          os << " \033[30;4m" << from_ << "\033[0m\n";
        os << "\t" << message_.str();
        return os.str();
      }

    private:
      inline static char* now() {
        static char buffer[10];
        time_t rawtime;
        time( &rawtime );
        struct tm* timeinfo = localtime( &rawtime );
        strftime( buffer, 10, "%H:%M:%S", timeinfo );
        return buffer;
      }
      /// Extract a full exception message
      inline std::string fullMessage() const {
        if ( type_ == Type::info || type_ == Type::debug || type_ == Type::warning )
          return shortMessage();
        if ( type_ == Type::verbatim )
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
      Type type_;
      /// Integer exception number
      int error_num_;
  };
  /// \brief Placeholder for debugging messages if logging threshold is not reached
  /// \date Apr 2018
  struct NullStream
  {
    /// Construct from a module name
    explicit NullStream( const char* ) {}
    /// Construct from a module name
    explicit NullStream( const std::string& ) {}
    /// Copy constructor
    NullStream( const Exception& ) {}
    /// Stream operator (null and void)
    template<class T> NullStream& operator<<( const T& ) { return *this; }
  };
}

#endif

