#ifndef CepGen_Modules_EventModifier_h
#define CepGen_Modules_EventModifier_h

#include <string>
#include <vector>

namespace cepgen
{
  class Event;
  class Parameters;
  class ParametersList;
  /// Class template to interface (external/internal) events modification algorithms
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date July 2019
  class EventModifier
  {
    public:
      /// Default constructor for an undefined modifier
        /// \param[in] params User-controlled steering parameters for this module
      explicit EventModifier( const ParametersList& params );
      /// Virtual destructor
      virtual ~EventModifier() = default;

      /// Set all runtime parameters steering this module
      virtual void setParameters( const Parameters& params ) { params_ = &params; }
      /// \brief Specify a random numbers generator seed for the external module
      /// \param[in] seed A RNG seed
      void setSeed( long long seed ) { seed_ = seed; }
      /// Return a human-readable name for this modifier
      const std::string& name() const { return name_; }

      /// Parse a configuration string
      virtual void readString( const char* ) {}
      /// Parse a configuration string
      virtual void readString( const std::string& param ) { readString( param.c_str() ); }
      /// Parse a list of configuration strings
      virtual void readStrings( const std::vector<std::string>& params );

      /// Initialise the event modifier before its running
      virtual void init() = 0;
      /** \brief Modify a full event
       * \param[inout] ev Input/output event
       * \param[inout] weight Event weight after modification
       * \param[in] full Perform the full state modification
       * \return Boolean stating whether or not the modification occured successfully
       */
      virtual bool run( Event& ev, double& weight, bool full ) = 0;
      /// Specify the process cross section, in pb
      virtual void setCrossSection( double xsec, double xsec_err ) {}

    protected:
      /// Name of the algorithm
      std::string name_;
      /// Random numbers generator seed fed to the algorithm
      long long seed_;
      /// Maximal number of trials for the algorithm
      unsigned short max_trials_;
      /// List of runtime parameters steering this module
      const Parameters* params_; // not owning
  };
}

#endif
