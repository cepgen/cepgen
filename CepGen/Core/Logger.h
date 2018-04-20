#ifndef CepGen_Core_Logger_h
#define CepGen_Core_Logger_h

#include <iostream>
#include <vector>
#include <regex>

namespace CepGen
{
  /**
   * \brief General purposes logger
   * \author Laurent Forthomme <laurent.forthomme@cern.ch>
   * \date 15 Oct 2015
   */
  class Logger
  {
    public:
      /// Logging threshold for the output stream
      enum class Level { Nothing = 0, Error, Warning, Information, Debug, DebugInsideLoop };

    private:
      /// Initialize a logging object
      Logger( std::ostream* os = &std::cout ) : level( Level::Information ), output( os ) {}
      ~Logger() {
        if ( output.get() == &std::cout )
          output.reset();
      }

      std::vector<std::regex> allowed_exc_;

    public:
      /// Retrieve the running instance of the logger
      static Logger& get( std::ostream* os = &std::cout ) {
        static Logger log( os );
        return log;
      }
      /// Add a new rule to display exceptions/messages
      /// \param[in] rule Regex rule to handle
      void addExceptionRule( const char* rule );
      /// Is the module set to be displayed/logged?
      /// \param[in] tmpl Module name to probe
      bool passExceptionRule( const std::string& tmpl ) const;

      /// Redirect the logger to a given output stream
      friend std::ostream& operator<<( std::ostream& os, const Logger::Level& lvl ) {
        switch ( lvl ) {
          case Level::Nothing:         return os << "None";
          case Level::Error:           return os << "Errors";
          case Level::Warning:         return os << "Warnings";
          case Level::Information:     return os << "Infos";
          case Level::Debug:           return os << "Debug";
          case Level::DebugInsideLoop: return os << "Debug (in loops)";
        }
        return os;
      }
      /// Logging threshold for the output stream
      Level level;
      /// Output stream to use for all logging operations
      std::shared_ptr<std::ostream> output;
  };
}

#endif
