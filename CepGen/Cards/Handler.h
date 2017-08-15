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
  /// Location for all steering card parsers/writers
  namespace Cards
  {
    /// List of all steering card types handled
    enum Type {
      Lpair, ///< LPAIR steering card
      Tcl
    };

    /// Generic steering card handler
    template<Type T>
    class Handler
    {
      public:
        /// Build a configuration from an external steering card
        /// \param[in] file Input file to parse
        Handler( const char* file );
        ~Handler() {}

        /// Store a configuration into an external steering card
        void store( const char* file ) {}
        /// Retrieve a configuration from a parsed steering cart
        Parameters& parameters() { return params_; }

      private:
        Parameters params_;
    };

    /// LPAIR steering cards handler
    typedef Handler<Lpair> LpairReader;
  }
}

#endif
