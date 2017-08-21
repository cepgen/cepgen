#ifndef CepGen_Cards_Handler_h
#define CepGen_Cards_Handler_h

#include "CepGen/Parameters.h"

#include "CepGen/Processes/GamGamLL.h"
#include "CepGen/Processes/PPtoLL.h"

#include "CepGen/Hadronisers/Pythia6Hadroniser.h"

namespace CepGen
{
  /// Location for all steering card parsers/writers
  namespace Cards
  {
    /// Generic steering card handler
    class Handler
    {
      public:
        /// Build a configuration from an external steering card
        Handler() {}
        ~Handler() {}

        /// Retrieve a configuration from a parsed steering cart
        Parameters& parameters() { return params_; }

      protected:
        Parameters params_;
    };
  }
}

#endif
