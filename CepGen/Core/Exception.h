#ifndef CepGen_Core_Exception_h
#define CepGen_Core_Exception_h

#include <csignal>

#include "CepGen/Utils/Logger.h"

namespace cepgen
{
  /// A generic exception type
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date 27 Mar 2015
  struct Exception
  {
    explicit inline Exception() = default;
    virtual ~Exception() noexcept = default;
    /// Enumeration of exception severities
    enum class Type {
      undefined = -1, ///< Irregular exception
      debug, ///< Debugging information to be enabled
      verbatim, ///< Raw information
      info, ///< Prettified information
      warning, ///< Casual non-stopping warning
      error, ///< General non-stopping error
      fatal ///< Critical and stopping error
    };
    /// Stream operator (null and void)
    template<class T> Exception& operator<<( const T& ) { return *this; }
    /// Dump the full exception information in a given output stream
    /// \param[inout] os the output stream where the information is dumped
    virtual void dump( std::ostream& os = *utils::Logger::get().output ) const = 0;
    /// Exception message
    virtual std::string message() const = 0;
  };

  /// A simple exception handler
  /// \date 24 Mar 2015
  class LoggedException : public Exception
  {
    public:
      /// Generic constructor
      /// \param[in] module exception classifier
      /// \param[in] type exception type
      /// \param[in] id exception code (useful for logging)
      explicit inline LoggedException( const char* module = "", Type type = Type::undefined, short id = 0 ) :
        module_( module ), type_( type ), error_num_( id ) {}
      /// Generic constructor
      /// \param[in] from method invoking the exception
      /// \param[in] module exception classifier
      /// \param[in] type exception type
      /// \param[in] id exception code (useful for logging)
      explicit inline LoggedException( const char* from, const char* module, Type type = Type::undefined, short id = 0 ) :
        from_( from ), module_( module ), type_( type ), error_num_( id ) {}
      /// Copy constructor
      inline LoggedException( const LoggedException& rhs ) :
        from_( rhs.from_ ), module_( rhs.module_ ), message_( rhs.message_.str() ),
        type_( rhs.type_ ), error_num_( rhs.error_num_ ) {}
      /// Default destructor (potentially killing the process)
      inline ~LoggedException() noexcept override {
        if ( type_ != Type::undefined )
          dump();
        // we stop this process' execution on fatal exception
        if ( type_ == Type::fatal && raise( SIGINT ) != 0 )
          exit( 0 );
      }

      //----- Overloaded stream operators

      /// Generic templated message feeder operator
      template<typename T>
      inline friend const LoggedException& operator<<( const LoggedException& exc, T var ) {
        LoggedException& nc_except = const_cast<LoggedException&>( exc );
        nc_except.message_ << var;
        return exc;
      }
      /// Generic templated vector-variables feeder operator
      template<typename T>
      inline friend const LoggedException& operator<<( const LoggedException& exc, std::vector<T> vec_var ) {
        LoggedException& nc_except = const_cast<LoggedException&>( exc );
        std::string sep( "{" );
        for ( const auto& var : vec_var )
          nc_except.message_ << sep << var, sep = ", ";
        nc_except.message_ << "}";
        return exc;
      }
      /// Pipe modifier operator
      inline friend const LoggedException& operator<<( const LoggedException& exc, std::ios_base&( *f )( std::ios_base& ) ) {
        LoggedException& nc_except = const_cast<LoggedException&>( exc );
        f( nc_except.message_ );
        return exc;
      }

      inline std::string message() const override {
        return message_.str();
      }

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
          case Type::debug: return "\033[33;1m Debug \033[0m";
          case Type::error: return "\033[31;1m Error \033[0m";
          case Type::fatal: return "\033[31;1m Fatal \033[0m";
          case Type::undefined: default: return "\33[7;1mUndefined\033[0m";
        }
      }

      inline void dump( std::ostream& os = *utils::Logger::get().output ) const override {
        if ( !utils::Logger::get().output )
          return;
        os << fullMessage() << std::endl;
      }
      /// Extract a one-line summary of the exception
      inline std::string shortMessage() const {
        std::ostringstream os;
        os << "[" << typeString() << "]";
        if ( type_ == Type::warning || type_ == Type::debug )
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
        switch ( type_ ) {
          case Type::info:
          case Type::debug:
          case Type::warning:
            return shortMessage();
          case Type::verbatim:
            return message_.str();
          default: {
            std::ostringstream os;
            os << "=============================== " << typeString() << " =============================== "
               << now() << std::endl;
            if ( !from_.empty() )
              os << " Raised by: " << from_ << "\n"
                 << " " << message_.str() << std::endl;
            if ( errorNumber() != 0 )
              os << "--------------------------------------------------------------------------------"
                 << "\n Error #" << error_num_ << std::endl;
            os << "================================================================================";
            return os.str();
          }
        }
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
      short error_num_;
  };

  /// Placeholder for debugging messages if logging threshold is not reached
  /// \date Apr 2018
  struct NullStream : Exception
  {
    using Exception::Exception;
    /// Empty constructor
    inline NullStream() {}
    /// Empty constructor
    inline NullStream( const LoggedException& ) {}
    void dump( std::ostream& os = *utils::Logger::get().output ) const override {}
    std::string message() const override { return ""; }
  };
}

#ifdef _WIN32
# define __FUNC__ __FUNCSIG__
#else
# define __FUNC__ __PRETTY_FUNCTION__
#endif

#define CG_LOG( mod ) \
  ( !CG_LOG_MATCH( mod, information ) ) \
  ? cepgen::NullStream() \
  : cepgen::LoggedException( __FUNC__, mod, cepgen::Exception::Type::verbatim )
#define CG_INFO( mod ) \
  ( !CG_LOG_MATCH( mod, information ) ) \
  ? cepgen::NullStream() \
  : cepgen::LoggedException( __FUNC__, mod, cepgen::Exception::Type::info )
#define CG_DEBUG( mod ) \
  ( !CG_LOG_MATCH( mod, debug ) ) \
  ? cepgen::NullStream() \
  : cepgen::LoggedException( __FUNC__, mod, cepgen::Exception::Type::debug )
#define CG_DEBUG_LOOP( mod ) \
  ( !CG_LOG_MATCH( mod, debugInsideLoop ) ) \
  ? cepgen::NullStream() \
  : cepgen::LoggedException( __FUNC__, mod, cepgen::Exception::Type::debug )
#define CG_WARNING( mod ) \
  ( !CG_LOG_MATCH( mod, warning ) ) \
  ? cepgen::NullStream() \
  : cepgen::LoggedException( __FUNC__, mod, cepgen::Exception::Type::warning )
#define CG_ERROR( mod ) \
  cepgen::LoggedException( __FUNC__, mod, cepgen::Exception::Type::error )
#define CG_FATAL( mod ) \
  cepgen::LoggedException( __FUNC__, mod, cepgen::Exception::Type::fatal )

#endif
