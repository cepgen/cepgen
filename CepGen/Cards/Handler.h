#ifndef CepGen_Cards_Handler_h
#define CepGen_Cards_Handler_h

#include <fstream>
#include <string>

#include "CepGen/Parameters.h"

#include "CepGen/Processes/GamGamLL.h"
#include "CepGen/Processes/PPtoLL.h"

#include "CepGen/Hadronisers/Pythia6Hadroniser.h"
#include "CepGen/Hadronisers/Jetset7Hadroniser.h"
#include "CepGen/Hadronisers/Herwig6Hadroniser.h"

namespace CepGen
{
  namespace Cards
  {
    enum Type { LPAIR };

    template<Type T>
    class Handler
    {
      public:
        Handler( const char* file );
        ~Handler();

        void parse( const char* file ) {}
        void store( const char* file ) {}
        Parameters& parameters() { return params_; }

      private:
        Parameters params_;
    };

    typedef Handler<LPAIR> LpairReader;
  }
}

#endif
