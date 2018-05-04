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
      Logger( std::ostream* os = &std::cout );
      ~Logger();

      std::vector<std::regex> allowed_exc_;

    public:
      /// Retrieve the running instance of the logger
      static Logger& get( std::ostream* os = &std::cout ) {
        static Logger log( os );
        return log;
      }
      /// Add a new rule to display exceptions/messages
      /// \param[in] rule Regex rule to handle
      void addExceptionRule( const std::string& rule );
      /// Is the module set to be displayed/logged?
      /// \param[in] tmpl Module name to probe
      bool passExceptionRule( const std::string& tmpl, const Level& lev ) const;

      /// Redirect the logger to a given output stream
      friend std::ostream& operator<<( std::ostream& os, const Logger::Level& lvl );
      /// Logging threshold for the output stream
      Level level;
      /// Output stream to use for all logging operations
      std::ostream* output;
  };

#if defined(__CINT__) || defined(__CLING__)
  inline Logger::Logger( std::ostream* ) {}
  inline Logger::~Logger() {}
  inline void Logger::addExceptionRule( const std::string& ) {}
  inline bool Logger::passExceptionRule( const std::string&, const Logger::Level& ) const { return false; }
#endif
}

#endif

