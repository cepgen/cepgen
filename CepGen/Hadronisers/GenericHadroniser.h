#ifndef CepGen_Hadronisers_GenericHadroniser_h
#define CepGen_Hadronisers_GenericHadroniser_h

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
     * Class template to define any hadroniser as a general object with defined methods
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date January 2014
     */
    class GenericHadroniser
    {
      public:
        friend std::ostream& operator<<( std::ostream& os, const GenericHadroniser& hadr ) { os << hadr.name().c_str(); return os; }
        friend std::ostream& operator<<( std::ostream& os, const GenericHadroniser* hadr ) { os << hadr->name().c_str(); return os; }

        /// Default constructor for an undefined hadroniser
        explicit GenericHadroniser( const char* name = "unnamed_hadroniser" ) : name_( name ) {}
        virtual ~GenericHadroniser() {}

        /// Hadronise a full event
        /// \param[inout] ev Event to hadronise
        /// \param[inout] weight Event weight after hadronisation
        /// \param[in] full Perform the full state hadronisation (incl. remnants fragmentation)
        /// \return Boolean stating whether or not the hadronisation occured successfully
        virtual bool run( Event& ev, double& weight, bool full ) = 0;
        /// Specify a random numbers generator seed for the hadroniser
        /// \param[in] seed A RNG seed
        virtual void setSeed( long long seed ) = 0;
        virtual void setCrossSection( double xsec, double xsec_err ) {};

        /// Return a human-readable name for this hadroniser
        inline std::string name() const { return name_; }

      protected:
        /// Name of the hadroniser
        std::string name_;
    };
  }
}

#endif
