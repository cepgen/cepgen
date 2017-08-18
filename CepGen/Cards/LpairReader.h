#ifndef CepGen_Cards_LpairReader_h
#define CepGen_Cards_LpairReader_h

#include "Handler.h"

#include <fstream>
#include <string>

namespace CepGen
{
  namespace Cards
  {
    class LpairReader : public Handler
    {
      public:
        LpairReader( const char* file );

        void store( const char* file ) const;

      private:
        //void parseKinematics( const libconfig::Setting& );
    };
  }
}

#endif
