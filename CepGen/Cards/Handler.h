#ifndef CepGen_Cards_Handler_h
#define CepGen_Cards_Handler_h

#include <fstream>
#include <string>

#include "CepGen/Parameters.h"

#include "CepGen/Processes/GamGamLL.h"
#include "CepGen/Processes/PPtoLL.h"

#include "CepGen/Hadronisers/Pythia6Hadroniser.h"

namespace CepGen
{
  namespace Cards
  {
    enum Type { Lpair, Tcl };

    template<Type T>
    class Handler
    {
      public:
        Handler( const char* file );
        ~Handler() {}

        void store( const char* file ) {}
        Parameters& parameters() { return params_; }

      private:
        Parameters params_;
    };

    typedef Handler<Lpair> LpairReader;
  }
}

#endif
