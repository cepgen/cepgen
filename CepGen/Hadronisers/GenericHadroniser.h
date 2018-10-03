#ifndef CepGen_Hadronisers_GenericHadroniser_h
#define CepGen_Hadronisers_GenericHadroniser_h

#include <vector>
#include <memory>
#include <iostream>

namespace CepGen
{
  class Event;
  class Particle;
  class Parameters;
  /// Location for all hadronisers to be run downstream to the events generation
  namespace Hadroniser
  {
    /**
     * \brief Class template to define any hadroniser as a general object with defined methods
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date January 2014
     */
    class GenericHadroniser
    {
      public:
        friend std::ostream& operator<<( std::ostream& os, const GenericHadroniser& hadr );
        friend std::ostream& operator<<( std::ostream& os, const GenericHadroniser* hadr );

        /// Default constructor for an undefined hadroniser
        explicit GenericHadroniser( const char* name = "unnamed_hadroniser" );
        virtual ~GenericHadroniser() {}

        /// Parse a configuration string
        virtual void readString( const char* ) {}
        /// Parse a configuration string
        virtual void readString( const std::string& param ) { readString( param.c_str() ); }
        /// Parse a list of configuration strings
        virtual void readStrings( const std::vector<std::string>& params );

        virtual void init() = 0;
        /// Hadronise a full event
        /// \param[inout] ev Event to hadronise
        /// \param[inout] weight Event weight after hadronisation
        /// \param[in] full Perform the full state hadronisation (incl. remnants fragmentation)
        /// \return Boolean stating whether or not the hadronisation occured successfully
        virtual bool run( Event& ev, double& weight, bool full ) = 0;
        /// Specify a random numbers generator seed for the hadroniser
        /// \param[in] seed A RNG seed
        virtual void setSeed( long long seed ) = 0;
        /// Specify the process cross section, in pb
        virtual void setCrossSection( double xsec, double xsec_err ) {}

        /// Return a human-readable name for this hadroniser
        std::string name() const;

      protected:
        /// Name of the hadroniser
        std::string name_;
    };
  }
}

#endif
