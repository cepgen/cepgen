#ifndef CepGen_Core_Logger_h
#define CepGen_Core_Logger_h

#include <iostream>

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
      enum LoggingLevel { Nothing=0, Error, Warning, Information, Debug, DebugInsideLoop };

    private:
      /// Initialize a logging object
      Logger() : level( Warning ), outputStream( std::cout ) {}
      ~Logger() {}

    public:
      /// Retrieve the running instance of the logger
      static Logger& get() {
        static Logger log;
        return log;
      }

      /// Redirect the logger to a given output stream
      friend std::ostream& operator<<( std::ostream& os, const Logger::LoggingLevel& lvl ) {
        switch ( lvl ) {
          case Logger::Nothing:         os << "None"; break;
          case Logger::Error:           os << "Errors"; break;
          case Logger::Warning:         os << "Warnings"; break;
          case Logger::Information:     os << "Infos"; break;
          case Logger::Debug:           os << "Debug"; break;
          case Logger::DebugInsideLoop: os << "Debug (in loops)"; break;
        }
        return os;
      }
      /// Logging threshold for the output stream
      LoggingLevel level;
      /// Output stream to use for all logging operations
      std::ostream& outputStream;
  };
}
  
#endif
