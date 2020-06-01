#ifndef CepGen_Cards_Handler_h
#define CepGen_Cards_Handler_h

#include "CepGen/Parameters.h"

namespace cepgen
{
  class Parameters;
  class ParametersList;
  /// Location for all steering card parsers/writers
  namespace card
  {
    /// Base steering card module
    class Handler
    {
      public:
        /// Build a configuration from an external steering card
        explicit Handler( const ParametersList& );
        ~Handler() = default;

        /// Small utility to retrieve the extension of a filename
        ///  (naive approach)
        static std::string extension( const std::string& file ) {
          return file.substr( file.find_last_of( "." )+1 );
        }
        /// Get the list of runtime parameters parsed
        Parameters* parameters() { return params_; }
        /// Specify runtime parameters to the handler
        virtual void pack( const Parameters* ) {};
        /// Retrieve a configuration from a parsed steering card
        virtual Parameters* parse( const std::string& filename, Parameters* params ) {
          return params;
        }
        /// Build a configuration from a steering card
        static Parameters* parse( const std::string& filename );
        /// Write the current configuration into a steering card
        virtual void write( const std::string& filename ) const {}
        /// Write a steering card from a configuration
        static void write( const Parameters* params, const std::string& filename );

      protected:
        /// Input filename
        const std::string filename_;
        /// List of parameters parsed from a card handler
        Parameters* params_;
    };
  }
}

#endif
