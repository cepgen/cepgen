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
    Logger() : Level(Warning), OutputStream(std::cout) {;}
    /// The static object present everywhere at runtime.
    static Logger* fLogger;
    /// A boolean stating whether or not the static object is already
    /// built.
    static bool fBuilt;

  public:    
    ~Logger() { fBuilt=false; }
    static Logger* GetInstance();
    
    enum LoggingLevel { Nothing, Error, Warning, Information, Debug, DebugInsideLoop };
    LoggingLevel Level;
    std::ostream& OutputStream;

  private:
};
  
#endif
