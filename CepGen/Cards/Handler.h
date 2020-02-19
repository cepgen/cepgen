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
        Handler() = default;
        ~Handler() = default;

        Parameters& parameters() { return params_; }

        /// Retrieve a configuration from a parsed steering card
        static std::unique_ptr<Handler> parse( const std::string& filename );

      protected:
        /// Small utility to retrieve the extension of a filename
        ///  (naive approach)
        static std::string extension( const std::string& file ) {
          return file.substr( file.find_last_of( "." )+1 );
        }
        static constexpr const char* FILENAME_KEY = "filename";
        /// List of parameters parsed from a card handler
        Parameters params_;
    };
  }
}

#endif
