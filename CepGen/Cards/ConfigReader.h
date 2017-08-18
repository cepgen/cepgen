#ifndef CepGen_Cards_ConfigReader_h
#define CepGen_Cards_ConfigReader_h

#include "Handler.h"

#ifdef LIBCONFIG
#include <libconfig.h++>
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

        void store( const char* file ) const;

      private:
        void parseIncomingKinematics( const libconfig::Setting& );
        void parseOutgoingKinematics( const libconfig::Setting& );
        void parseVegas( const libconfig::Setting& );
        void parseGenerator( const libconfig::Setting& );
#else
      public:
        inline ConfigReader( const char* ) { InWarning( "libconfig++ is not present on this machine" ); }
#endif
    };
  }
}

#endif
