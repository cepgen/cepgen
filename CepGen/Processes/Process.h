#ifndef CepGen_Processes_Process_h
#define CepGen_Processes_Process_h

#include "CepGen/Core/ParametersList.h"

#include "CepGen/Physics/Kinematics.h"
#include "CepGen/Event/Particle.h"

#include <map>
#include <vector>
#include <memory>

namespace cepgen
{
  class FormFactors;
  class Event;
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
        /// \param[in] has_event Do we generate the associated event structure?
        explicit Process( const ParametersList& params, bool has_event = true );
        /// Copy constructor for a user process
        Process( const Process& );
        virtual ~Process() = default;

        /// Reset process prior to the phase space and variables definition
        void clear();

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

      public:
        /// Copy all process attributes into a new object
        virtual std::unique_ptr<Process> clone( const ParametersList& params = ParametersList() ) const;
        /// Set the incoming and outgoing state to be expected in the process
        inline virtual void addEventContent() {}
        /// Compute the phase space point weight
        virtual double computeWeight() = 0;
        /// Compute the incoming state kinematics
        virtual void prepareKinematics() = 0;
        /// Fill the Event object with the particles' kinematics
        /// \param[in] symmetrise Symmetrise the event? (randomise the production of positively- and negatively-charged outgoing central particles)
        virtual void fillKinematics( bool symmetrise = false ) = 0;

      public:
        /// Restore the Event object to its initial state
        void clearEvent();
        /// Set the list of kinematic cuts to apply on the outgoing particles' final state
        /// \param[in] kin The Kinematics object containing the kinematic parameters
        void setKinematics( const Kinematics& kin );
        /**
         * Sets the phase space point to compute the weight associated to it.
         * \brief Sets the phase space point to compute
         * \param[in] ndim The number of dimensions of the point in the phase space
         * \param[in] x The (\a ndim_)-dimensional point in the phase space on which the kinematics and the cross-section are computed
         */
        void setPoint( double* x, const size_t ndim );
        /// Compute the weight for this point in the phase-space
        double weight();
        /// Dump the evaluated point's coordinates in the standard output stream
        void dumpPoint() const;
        /// List all variables handled by this generic process
        void dumpVariables() const;

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
        bool hasEvent() const { return (bool)event_; }
        /// Event object containing all the generated Particle objects and their relationships
        inline const Event& event() const { return *event_; }
        /// Non-const event retrieval method
        inline Event& event() { return *event_; }
        /// Event pointer retrieval method
        inline Event* eventPtr() { return event_.get(); }

      protected:
        const double mp_; ///< Proton mass, in GeV/c\f$^2\f$
        const double mp2_; ///< Squared proton mass, in GeV\f$^2\f$/c\f$^4\f$
        /// Type of mapping to apply on the variable
        enum class Mapping
        {
          linear = 0, ///< a linear \f${\rm d}x\f$ mapping
          exponential, ///< an exponential \f$\frac{\dot{x}}{x} = \dot{\log x}\f$ mapping
          square, ///< a square \f${\rm d}x^2=2x\cdot\dot{x}\f$ mapping
          power_law ///< a power-law mapping inherited from LPAIR
          /**
           * Define modified variables of integration to avoid peaks integrations (see \cite Vermaseren:1982cz for details):
           * - \f$y_{\rm out} = x_{\rm min}\left(\frac{x_{\rm max}}{x_{\rm min}}\right)^{\rm exp}\f$ the new variable
           * - \f${\rm d}y_{\rm out} = x_{\rm min}\left(\frac{x_{\rm max}}{x_{\rm min}}\right)^{\rm exp}\log\frac{x_{\rm min}}{x_{\rm max}}\f$, the new variable's differential form
           * \note This method overrides the set of `mapxx` subroutines in ILPAIR, with a slight difference according to the sign of the
           *  \f${\rm d}y_{\rm out}\f$ parameter :
           *  - left unchanged :
           * > `mapw2`, `mapxq`, `mapwx`, `maps2`
           *  - opposite sign :
           * > `mapt1`, `mapt2`
           */
        };
        /// Human-friendly printout of the mapping type
        friend std::ostream& operator<<( std::ostream&, const Mapping& );
        /// Register a variable to be handled and populated whenever
        ///  a new phase space point weight is to be calculated.
        /// \note To be run once per generation (before any point computation)
        /// \param[out] out Reference to the variable to be mapped
        /// \param[in] type Type of mapping to apply
        /// \param[in] in Integration limits
        /// \param[in] default_limits Limits to apply if none retrieved from the user configuration
        /// \param[in] description Human-readable description of the variable
        Process& defineVariable( double& out, const Mapping& type, Limits in = { 0., 1. }, const Limits& default_limits = { 0., 1. }, const std::string& description = "" );
        /// Generate and initialise all variables handled by this process
        /// \note To be run at each point computation (therefore, to be optimised!)
        void generateVariables() const;
        /// Phase space point-dependent component of the Jacobian weight of the point in the phase space for integration
        double jacobian() const;

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
        /// Numerical limits for sanity comparisons
        static constexpr double NUM_LIMITS = 1.e-3; // MeV/mm-level
        /// Handler to a variable mapped by this process
        struct MappingVariable
        {
          std::string description; ///< Human-readable description of the variable
          Limits limits; ///< Kinematic limits to apply on the variable
          double& value; ///< Reference to the process variable to generate/map
          Mapping type; ///< Interpolation type
          unsigned short index; ///< Corresponding integration variable
        };
        /// Collection of variables to be mapped at the weight generation stage
        std::vector<MappingVariable> mapped_variables_;
        /// Point coordinate for matrix element computation
        std::vector<double> point_coord_;
        /// Phase space point-independant component of the Jacobian weight of the point in the phase space for integration
        double base_jacobian_;
        /// \f$s\f$, squared centre of mass energy of the incoming particles' system, in \f$\mathrm{GeV}^2\f$
        double s_;
        /// \f$\sqrt s\f$, centre of mass energy of the incoming particles' system (in GeV)
        double sqs_;
        /// first incoming beam particle squared mass
        double mA2_;
        /// second incoming beam particle squared mass
        double mB2_;
        /// First diffractive state squared mass
        double mX2_;
        /// Second diffractive state squared mass
        double mY2_;
        /// First parton virtuality
        double t1_;
        /// Second parton virtuality
        double t2_;

        /// Set of cuts to apply on the final phase space
        Kinematics kin_;
        /// Event object containing all the information on all particles in the system
        std::unique_ptr<Event> event_;
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
