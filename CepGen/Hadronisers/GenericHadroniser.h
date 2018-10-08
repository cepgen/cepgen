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
  class ParametersList;
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
        /// Write out all hadroniser attributes in output stream
        friend std::ostream& operator<<( std::ostream& os, const GenericHadroniser& hadr );
        /// Write out all hadroniser attributes in output stream
        friend std::ostream& operator<<( std::ostream& os, const GenericHadroniser* hadr );

        /// Default constructor for an undefined hadroniser
        explicit GenericHadroniser( const char* name, const ParametersList& );
        virtual ~GenericHadroniser() {}

        /// Parse a configuration string
        virtual void readString( const char* ) {}
        /// Parse a configuration string
        virtual void readString( const std::string& param ) { readString( param.c_str() ); }
        /// Parse a list of configuration strings
        virtual void readStrings( const std::vector<std::string>& params );
        /// Initialise the event hadroniser before its running
        virtual void init() = 0;
        /** \brief Hadronise a full event
         * \param[inout] ev Event to hadronise
         * \param[inout] weight Event weight after hadronisation
         * \param[in] full Perform the full state hadronisation (incl. remnants fragmentation)
         * \return Boolean stating whether or not the hadronisation occured successfully
         */
        virtual bool run( Event& ev, double& weight, bool full ) = 0;
        /// Specify the process cross section, in pb
        virtual void setCrossSection( double xsec, double xsec_err ) {}

        /// \brief Specify a random numbers generator seed for the hadroniser
        /// \param[in] seed A RNG seed
        void setSeed( long long seed ) { seed_ = seed; }

        /// Return a human-readable name for this hadroniser
        std::string name() const;

      protected:
        /// Name of the hadroniser
        std::string name_;
        /// Random numbers generator seed for the hadroniser
        long long seed_;
        /// Maximal number of trials for the hadronisation of the proton(s) remnants
        unsigned short max_trials_;
    };
  }
}

#endif
