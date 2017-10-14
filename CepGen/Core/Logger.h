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
      static Logger& get();

      /// Redirect the logger to a given output stream
      friend std::ostream& operator<<( std::ostream& os, const Logger::LoggingLevel& lvl );
      /// Logging threshold for the output stream
      LoggingLevel level;
      /// Output stream to use for all logging operations
      std::ostream& outputStream;
  };
}
  
#endif
