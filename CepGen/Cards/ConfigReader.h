#ifndef CepGen_Cards_ConfigReader_h
#define CepGen_Cards_ConfigReader_h

#include "Handler.h"

#ifdef LIBCONFIG
#include <libconfig.h++>
#else
#pragma message "libconfig++ is not found on your system!"
#endif

namespace CepGen
{
  namespace Cards
  {
    class ConfigReader : public Handler
    {
#ifdef LIBCONFIG
      public:
        ConfigReader( const char* file );

        static void store( const Parameters*, const char* file );

      private:
        void parseIncomingKinematics( const libconfig::Setting& );
        void parseOutgoingKinematics( const libconfig::Setting& );
        void parseVegas( const libconfig::Setting& );
        void parseGenerator( const libconfig::Setting& );

        static void writeProcess( const Parameters*, libconfig::Setting& );
        static void writeIncomingKinematics( const Parameters*, libconfig::Setting& );
        static void writeOutgoingKinematics( const Parameters*, libconfig::Setting& );
        static void writeVegas( const Parameters*, libconfig::Setting& );
        static void writeGenerator( const Parameters*, libconfig::Setting& );
#else
      public:
        inline ConfigReader( const char* ) { InWarning( "libconfig++ is not present on this machine" ); }
        static void store( const Parameters*, const char* file ) {}
#endif
    };
  }
}

#endif
