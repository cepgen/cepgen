#ifndef CepGen_Cards_ConfigHandler_h
#define CepGen_Cards_ConfigHandler_h

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
    /// CepGen configuration cards reader/writer
    class ConfigHandler : public Handler
    {
      public:
        /// Read a standard configuration card
        ConfigHandler( const char* file );

        /// Store a configuration into a steering card
        static void store( const Parameters*, const char* file );

#ifdef LIBCONFIG
      private:
        void parseIncomingKinematics( const libconfig::Setting& );
        void parseOutgoingKinematics( const libconfig::Setting& );
        void parseIntegrator( const libconfig::Setting& );
        void parseGenerator( const libconfig::Setting& );
        void parseTamingFunctions( const libconfig::Setting& );
        void parseHadroniser( const libconfig::Setting& );

        static void writeProcess( const Parameters*, libconfig::Setting& );
        static void writeIncomingKinematics( const Parameters*, libconfig::Setting& );
        static void writeOutgoingKinematics( const Parameters*, libconfig::Setting& );
        static void writeTamingFunctions( const Parameters*, libconfig::Setting& );
        static void writeIntegrator( const Parameters*, libconfig::Setting& );
        static void writeHadroniser( const Parameters*, libconfig::Setting& );
        static void writeVegas( const Parameters*, libconfig::Setting& );
        static void writeGenerator( const Parameters*, libconfig::Setting& );
#endif
    };
  }
}

#endif
