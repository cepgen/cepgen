#ifndef CepGen_Processes_GenericProcess_h
#define CepGen_Processes_GenericProcess_h

#include "CepGen/Event/Event.h"
#include "CepGen/Physics/Kinematics.h"

#include <vector>
#include <memory>

namespace CepGen
{
  class FormFactors;
  /// Location for all physics processes to be generated
  namespace Process
  {
    /// Class template to define any process to compute using this MC integrator/events generator
    /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
    /// \date Jan 2014
    class GenericProcess
    {
      public:
        /// Default constructor for an undefined process
        /// \param[in] name Process name
        /// \param[in] description Human-readable description of the process
        /// \param[in] has_event Do we generate the associated event structure?
        GenericProcess( const std::string& name, const std::string& description = "<invalid process>", bool has_event = true );
        /// Copy constructor for a user process
        GenericProcess( const GenericProcess& );
        virtual ~GenericProcess() = default;

        /// Assignment operator
        GenericProcess& operator=( const GenericProcess& );

        /// Human-readable format dump of a GenericProcess object
        friend std::ostream& operator<<( std::ostream& os, const GenericProcess& proc );
        /// Human-readable format dump of a pointer to a GenericProcess object
        friend std::ostream& operator<<( std::ostream& os, const GenericProcess* proc );

        /// Generic map of particles with their role in the process
        typedef std::map<Particle::Role,PDG> ParticlesRoleMap;
        /// Pair of particle with their associated role in the process
        typedef std::pair<Particle::Role,PDG> ParticleWithRole;
        /// Map of all incoming state particles in the process
        typedef ParticlesRoleMap IncomingState;
        /// Map of all outgoing particles in the process
        typedef std::map<Particle::Role,std::vector<PDG> > OutgoingState;

        /// Copy all process' attributes into a new object
        virtual std::unique_ptr<GenericProcess> clone() const = 0;

        /// Restore the Event object to its initial state
        inline void clearEvent() { event_->restore(); }
        /// Set the kinematics of the incoming state particles
        void setIncomingKinematics( const Particle::Momentum& p1, const Particle::Momentum& p2 );
        /// Compute the incoming state kinematics
        void prepareKinematics();

      public:
        /// Set the incoming and outgoing state to be expected in the process
        inline virtual void addEventContent() {}
        /// Set the list of kinematic cuts to apply on the outgoing particles' final state
        /// \param[in] cuts The Cuts object containing the kinematic parameters
        virtual void setKinematics( const Kinematics& cuts );
        /// Return the number of dimensions on which the integration has to be performed
        /// \return Number of dimensions on which to integrate
        virtual unsigned int numDimensions( const Kinematics::Mode& ) const = 0;

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
        /// Complete list of Particle with their role in the process for the point considered in the phase space, returned as an Event object.
        /// \return Event object containing all the generated Particle objects
        inline std::shared_ptr<Event> event() const { return event_; }

        ///Get the number of dimensions on which the integration is performed
        inline const unsigned int ndim() const { return x_.size(); }
        /// Get the value of a component of the d-dimensional point considered
        double x( unsigned int idx ) const;
        /// Name of the process considered
        inline const std::string& name() const { return name_; }
        /// Human-readable description of the process
        inline const std::string& description() const { return description_; }

        /// Does the process contain (and hold) an event?
        bool hasEvent() const { return has_event_; }
        /// Pointer to the last event produced in this run
        std::shared_ptr<Event> last_event;

      protected:
        static const double mp_, mp2_;

        /// Set the incoming and outgoing states to be defined in this process (and prepare the Event object accordingly)
        void setEventContent( const IncomingState& ini, const OutgoingState& fin );
        /// Compute the electric/magnetic form factors for the two considered \f$Q^{2}\f$ momenta transfers
        void formFactors( double q1, double q2, FormFactors& fp1, FormFactors& fp2 ) const;

        // ---

        /// Name of the process
        std::string name_;
        /// Process human-readable description
        std::string description_;

      public:
        /// Is it the first time the process is computed?
        bool first_run;

      protected:
        /// Array of double precision floats representing the point on which the weight in the cross-section is computed
        std::vector<double> x_;
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
        Kinematics cuts_;
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
  }
  /// Helper typedef for a Process unique pointer
  typedef std::unique_ptr<Process::GenericProcess> ProcessPtr;
}

#endif
