#ifndef CepGen_Hadronisers_GenericHadroniser_h
#define CepGen_Hadronisers_GenericHadroniser_h

#include <memory>
#include <iostream>

namespace CepGen
{
  class Event;
  class Particle;
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
        GenericHadroniser( const char* name="unnamed_hadroniser" ) : name_( name ) {}
        virtual ~GenericHadroniser() {}

        /// Main caller to hadronise a single particle
        /// \param[in] part The Particle object which will be hadronised
        /// \return A boolean stating whether or not the hadronisation occured successfully
        virtual bool hadronise( const Particle& part, Event& ev ) = 0;
        /// Hadronise a full event
        /// \param[inout] ev Event to hadronise
        /// \return Boolean stating whether or not the hadronisation occured successfully
        virtual bool hadronise( Event& ev ) = 0;

        /// Return a human-readable name for this hadroniser
        inline std::string name() const { return name_; }

      protected:
        /// Name of the hadroniser
        std::string name_;
    };
  }
}

#endif
