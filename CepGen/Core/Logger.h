#ifndef CepGen_Core_Logger_h
#define CepGen_Core_Logger_h

#include <iostream>
#include <vector>
#include <regex>

namespace CepGen
{
  /// General purposes logger
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date 15 Oct 2015
  class Logger
  {
    public:
      /// Logging threshold for the output stream
      enum class Level { nothing = 0, error, warning, information, debug, debugInsideLoop };

    private:
      /// Initialize a logging object
      Logger( std::ostream* os = &std::cout ) :
        level( Level::information ), output( os ) {}

#if !defined(__CINT__) && !defined(__CLING__)
      std::vector<std::regex> allowed_exc_;
#endif

    public:
      /// Retrieve the running instance of the logger
      static Logger& get( std::ostream* os = &std::cout ) {
        static Logger log( os );
        return log;
      }
      /// Add a new rule to display exceptions/messages
      /// \param[in] rule Regex rule to handle
      void addExceptionRule( const std::string& rule ) {
#if !defined(__CINT__) && !defined(__CLING__)
        allowed_exc_.emplace_back( rule );
#endif
      }
      /// Is the module set to be displayed/logged?
      /// \param[in] tmpl Module name to probe
      bool passExceptionRule( const std::string& tmpl, const Level& lev ) const {
#if !defined(__CINT__) && !defined(__CLING__)
        if ( level >= lev )
          return true;
        if ( allowed_exc_.size() == 0 )
          return false;
        for ( const auto& rule : allowed_exc_ )
          if ( std::regex_match( tmpl, rule ) )
            return true;
#endif
        return false;
      }

      /// Redirect the logger to a given output stream
      friend std::ostream& operator<<( std::ostream& os, const Logger::Level& lvl ) {
        switch ( lvl ) {
          case Logger::Level::nothing:
            return os << "None";
          case Logger::Level::error:
            return os << "Errors";
          case Logger::Level::warning:
            return os << "Warnings";
          case Logger::Level::information:
            return os << "Infos";
          case Logger::Level::debug:
            return os << "Debug";
          case Logger::Level::debugInsideLoop:
            return os << "Debug (in loops)";
        }
        return os;
      }
      /// Logging threshold for the output stream
      Level level;
      /// Output stream to use for all logging operations
      std::ostream* output;
  };
}

#endif

