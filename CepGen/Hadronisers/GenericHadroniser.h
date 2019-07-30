#ifndef CepGen_Hadronisers_GenericHadroniser_h
#define CepGen_Hadronisers_GenericHadroniser_h

#include "CepGen/Core/EventModifier.h"

#include <vector>
#include <iosfwd>

namespace cepgen
{
  class Event;
  class Particle;
  class Parameters;
  class ParametersList;
  /// Location for all hadronisers to be run downstream to the events generation
  namespace hadr
  {
    /**
     * \brief Class template to define any hadroniser as a general object with defined methods
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date January 2014
     */
    class GenericHadroniser : public EventModifier
    {
      public:
        /// Write out all hadroniser attributes in output stream
        friend std::ostream& operator<<( std::ostream& os, const GenericHadroniser& hadr );
        /// Write out all hadroniser attributes in output stream
        friend std::ostream& operator<<( std::ostream& os, const GenericHadroniser* hadr );

        /// Default constructor for an undefined hadroniser
        explicit GenericHadroniser( const ParametersList&, const std::string& name = "<invalid hadroniser>" );

        /// Parse a configuration string
        virtual void readString( const char* ) {}
        /// Parse a configuration string
        virtual void readString( const std::string& param ) { readString( param.c_str() ); }
        /// Parse a list of configuration strings
        virtual void readStrings( const std::vector<std::string>& params );

        /// Specify whether the beam remnants are to be fragmented
        bool fragmentRemnants() const { return remn_fragm_; }

      protected:
        /// Switch on/off the remnants fragmentation where applicable
        const bool remn_fragm_;
    };
  }
}

#endif
