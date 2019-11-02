#ifndef CepGen_Core_Process_h
#define CepGen_Core_Process_h

#include "CepGen/Event/Particle.h"
#include "CepGen/Physics/Kinematics.h"
#include "CepGen/Core/ParametersList.h"

#include <map>
#include <vector>
#include <memory>

namespace cepgen
{
  class Event;
  class FormFactors;
  /// Location for all physics processes to be generated
  namespace proc
  {
    /// \brief Class template to define any process to compute using this MC integrator/events generator
    /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
    /// \date Jan 2014
    class Process
    {
      public:
        /// Default constructor for an undefined process
        /// \param[in] params Process-level parameters
        /// \param[in] name Process name
        /// \param[in] description Human-readable description of the process
        /// \param[in] has_event Do we generate the associated event structure?
        Process( const ParametersList& params, const std::string& name = "<invalid name>", const std::string& description = "<invalid process>", bool has_event = true );
        /// Copy constructor for a user process
        Process( const Process& );
        virtual ~Process() {}

        /// Assignment operator
        Process& operator=( const Process& );

        /// Human-readable format dump of a Process object
        friend std::ostream& operator<<( std::ostream& os, const Process& proc );
        /// Human-readable format dump of a pointer to a Process object
        friend std::ostream& operator<<( std::ostream& os, const Process* proc );

        /// Map of all incoming state particles in the process
        typedef std::map<Particle::Role,pdgid_t> IncomingState;
        /// Map of all outgoing particles in the process
        typedef std::map<Particle::Role,std::vector<pdgid_t> > OutgoingState;

        /// Copy all process attributes into a new object
        virtual std::unique_ptr<Process> clone( const ParametersList& params = ParametersList() ) const = 0;

        /// Restore the Event object to its initial state
        void clearEvent();
        /// Set the kinematics of the incoming state particles
        void setIncomingKinematics( const Momentum& p1, const Momentum& p2 );
        /// Compute the incoming state kinematics
        virtual void prepareKinematics() {}

      public:
        /// Set the incoming and outgoing state to be expected in the process
        inline virtual void addEventContent() {}
        /// Set the list of kinematic cuts to apply on the outgoing particles' final state
        /// \param[in] kin The Kinematics object containing the kinematic parameters
        virtual void setKinematics( const Kinematics& kin );
        /// Return the number of dimensions on which the integration has to be performed
        /// \return Number of dimensions on which to integrate
        virtual unsigned int numDimensions() const { return mapped_variables_.size(); }

        /// Prepare the process for its integration over the whole phase space
        inline virtual void beforeComputeWeight() {}
        /// Compute the weight for this point in the phase-space
        virtual double computeWeight() = 0;
        /// Fill the Event object with the particles' kinematics
        /// \param[in] symmetrise Symmetrise the event? (randomise the production of positively- and negatively-charged outgoing central particles)
        virtual void fillKinematics( bool symmetrise = false ) = 0;

      public:
        /**
         * Sets the phase space point to compute the weight associated to it.
         * \brief Sets the phase space point to compute
         * \param[in] ndim The number of dimensions of the point in the phase space
         * \param[in] x[] The (\a ndim_)-dimensional point in the phase space on which the kinematics and the cross-section are computed
         */
        void setPoint( const unsigned int ndim, double* x );
        /// Dump the evaluated point's coordinates in the standard output stream
        void dumpPoint() const;
        /// List all variables handled by this generic process
        void dumpVariables() const;
        /// Complete list of Particle with their role in the process for the point considered in the phase space, returned as an Event object.
        /// \return Event object containing all the generated Particle objects
        inline std::shared_ptr<Event> event() const { return event_; }

        ///Get the number of dimensions on which the integration is performed
        inline size_t ndim() const { return mapped_variables_.size(); }
        /// Get the value of a component of the d-dimensional point considered
        double x( unsigned int idx ) const;
        /// Process-specific parameters
        inline const ParametersList& parameters() const { return params_; }
        /// Name of the process considered
        inline const std::string& name() const { return name_; }
        /// Human-readable description of the process
        inline const std::string& description() const { return description_; }

        /// Does the process contain (and hold) an event?
        bool hasEvent() const { return has_event_; }
        /// Pointer to the last event produced in this run
        std::shared_ptr<Event> last_event;

      protected:
        const double mp_; ///< Proton mass, in GeV/c\f$^2\f$
        const double mp2_; ///< Squared proton mass, in GeV\f$^2\f$/c\f$^4\f$
        /// Type of mapping to apply on the variable
        enum class Mapping
        {
          /// a linear \f${\rm d}x\f$ mapping
          linear = 0,
          /// a logarithmic \f$\frac{{\rm d}x}{x} = {\rm d}(\log x)\f$ mapping
          logarithmic,
          /// a square \f${\rm d}x^2=2x\cdot{\rm d}x\f$ mapping
          square
        };
        friend std::ostream& operator<<( std::ostream&, const Mapping& );
        /// Register a variable to be handled and populated whenever
        ///  a new phase space point weight is to be calculated.
        /// \note To be run once per generation (before any point computation)
        /// \param[out] out Reference to the variable to be mapped
        /// \param[in] type Type of mapping to apply
        /// \param[in] in Integration limits
        /// \param[in] default_limits Limits to apply if none retrieved from the user configuration
        /// \param[in] description Human-readable description of the variable
        Process& defineVariable( double& out, const Mapping& type, const Limits& in, Limits default_limits, const char* description );
        /// Generate and initialise all variables handled by this process
        /// \return Phase space point-dependent component of the Jacobian weight of the point in the phase space for integration
        /// \note To be run at each point computation (therefore, to be optimised!)
        double generateVariables() const;

        /// Set the incoming and outgoing states to be defined in this process (and prepare the Event object accordingly)
        void setEventContent( const IncomingState& ini, const OutgoingState& fin );

        // ---

        /// Process-specific parameters
        ParametersList params_;
        /// Name of the process
        std::string name_;
        /// Process human-readable description
        std::string description_;

      public:
        /// Is it the first time the process is computed?
        bool first_run;

      protected:
        /// Handler to a variable mapped by this process
        struct MappingVariable
        {
          /// Human-readable description of the variable
          std::string description;
          /// Kinematic limits to apply on the variable
          Limits limits;
          /// Reference to the process variable to generate/map
          double& value;
          /// Interpolation type
          Mapping type;
          /// Corresponding integration variable
          unsigned short index;
        };
        /// Collection of variables to be mapped at the weight generation stage
        std::vector<MappingVariable> mapped_variables_;
        /// Phase space point-independant component of the Jacobian weight of the point in the phase space for integration
        double base_jacobian_;
        /// \f$s\f$, squared centre of mass energy of the incoming particles' system, in \f$\mathrm{GeV}^2\f$
        double s_;
        /// \f$\sqrt s\f$, centre of mass energy of the incoming particles' system (in GeV)
        double sqs_;
        /// Invariant mass of the first proton-like outgoing particle (or remnant)
        double MX_;
        /// Invariant mass of the second proton-like outgoing particle (or remnant)
        double MY_;
        /// \f$m_1^2\f$, squared mass of the first proton-like incoming particle
        double w1_;
        /// \f$m_2^2\f$, squared mass of the second proton-like incoming particle
        double w2_;
        /// Virtuality of the first incoming photon
        double t1_;
        /// Virtuality of the second incoming photon
        double t2_;

        /// Set of cuts to apply on the final phase space
        Kinematics kin_;
        /// Does the process contain (and hold) an event?
        bool has_event_;
        /// Event object containing all the information on the in- and outgoing particles
        std::shared_ptr<Event> event_;
        /// Is the phase space point set?
        bool is_point_set_;

      private:
        /**
         * Is the system's kinematics well defined and compatible with the process ?
         * This check is mandatory to perform the d-dimensional point's cross-section computation.
         * \brief Is the system's kinematics well defined?
         * \return A boolean stating if the input kinematics and the final states are well-defined
         */
        bool isKinematicsDefined();
    };
    /// Helper typedef for a Process unique pointer
    typedef std::unique_ptr<Process> ProcessPtr;
  }
}

#endif
