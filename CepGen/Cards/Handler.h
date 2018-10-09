#ifndef CepGen_Cards_Handler_h
#define CepGen_Cards_Handler_h

#include "CepGen/Parameters.h"

namespace cepgen
{
  class Parameters;
  /// Location for all steering card parsers/writers
  namespace card
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
        /// Small utility to retrieve the extension of a filename
        ///  (naive approach)
        static std::string getExtension( const char* filename ) {
          const std::string file( filename );
          return file.substr( file.find_last_of( "." )+1 );
        }

      protected:
        /// List of parameters parsed from a card handler
        Parameters params_;
    };
  }
}

#endif
