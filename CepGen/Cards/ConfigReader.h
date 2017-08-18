#ifndef CepGen_Cards_ConfigReader_h
#define CepGen_Cards_ConfigReader_h

#include "Handler.h"

#include <libconfig.h++>

namespace CepGen
{
  namespace Cards
  {
    class ConfigReader : public Handler
    {
      public:
        ConfigReader( const char* file );

        void store( const char* file ) const;

      private:
        void parseIncomingKinematics( const libconfig::Setting& );
        void parseOutgoingKinematics( const libconfig::Setting& );
        void parseVegas( const libconfig::Setting& );
        void parseGenerator( const libconfig::Setting& );
    };
  }
}

#endif
