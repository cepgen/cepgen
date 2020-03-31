#ifndef CepGen_Cards_Handler_h
#define CepGen_Cards_Handler_h

#include "CepGen/Parameters.h"

namespace cepgen
{
  class Parameters;
  /// Location for all steering card parsers/writers
  namespace card
  {
    /// Base steering card module
    class Handler
    {
      public:
        /// Build a configuration from an external steering card
        Handler() = default;
        ~Handler() = default;

        /// Small utility to retrieve the extension of a filename
        ///  (naive approach)
        static std::string extension( const std::string& file ) {
          return file.substr( file.find_last_of( "." )+1 );
        }
        /// Get the list of runtime parameters parsed
        Parameters& parameters() { return params_; }
        /// Retrieve a configuration from a parsed steering card
        virtual Parameters& parse( const std::string& filename, Parameters& params ) = 0;
        static Parameters& parse( const std::string& filename );

      protected:
        static constexpr const char* FILENAME_KEY = "filename";
        /// List of parameters parsed from a card handler
        Parameters params_;
    };
  }
}

#endif
