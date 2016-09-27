#ifndef Logger_h
#define Logger_h

#include <iostream>

/**
 * \brief General purposes logger
 * \author Laurent Forthomme <laurent.forthomme@cern.ch>
 * \date 15 Oct 2015 
 */
class Logger
{
  private:
    /// Initialize a logging object
    Logger() : Level( Warning ), OutputStream( std::cout ) {;}
    /// The static object present everywhere at runtime.
    static Logger* fLogger;
    /// A boolean stating whether or not the static object is already
    /// built.
    static bool fBuilt;

  public:
    ~Logger() { fBuilt=false; }
    /// Retrieve the running instance of the logger
    static Logger* GetInstance();
    
    /// Logging threshold for the output stream
    enum LoggingLevel { Nothing=0, Error, Warning, Information, Debug, DebugInsideLoop };
    /// Redirect the logger to a given output stream
    friend std::ostream& operator<<( std::ostream& os, const Logger::LoggingLevel& lvl );
    /// Logging threshold for the output stream
    LoggingLevel Level;
    /// Output stream to use for all logging operations
    std::ostream& OutputStream;

  private:
};
  
#endif
