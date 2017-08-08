#ifndef GenericProcess_h
#define GenericProcess_h

#include "CepGen/Physics/Event.h"
#include "CepGen/Physics/Kinematics.h"
#include "CepGen/Physics/Physics.h"
#include "CepGen/Physics/StructureFunctions.h"
#include "CepGen/Physics/FormFactors.h"

namespace CepGen
{
  namespace Process
  {
    /**
     * Class template to define any process to compute using this MC integrator/events generator
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Jan 2014
     */
    class GenericProcess
    {
      public:
        /// Human-readable format dump of a GenericProcess object
        friend std::ostream& operator<<( std::ostream& os, const GenericProcess& proc );
        /// Human-readable format dump of a pointer to a GenericProcess object
        friend std::ostream& operator<<( std::ostream& os, const GenericProcess* proc );

        /// Generic map of particles with their role in the process
        typedef std::map<Particle::Role,Particle::ParticleCode> ParticlesRoleMap;
        /// Pair of particle with their associated role in the process
        typedef std::pair<Particle::Role,Particle::ParticleCode> ParticleWithRole;
        /// Map of all incoming state particles in the process
        typedef ParticlesRoleMap IncomingState;
        /// Map of all outgoing particles in the process
        typedef ParticlesRoleMap OutgoingState;

        /// Default constructor for an undefined process
        /// \param[in] name Human-readable format of the process name
        GenericProcess( const std::string& name="<invalid process>" );
        virtual ~GenericProcess();

        /// Restore the Event object to its initial state
        inline void clearEvent() { event_->restore(); }
        /// Set the kinematics of the incoming state particles
        void setIncomingKinematics( const Particle::Momentum& p1, const Particle::Momentum& p2 );
        /// Compute the incoming state kinematics
        void prepareKinematics();

      public:
        /// Set the incoming and outgoing state to be expected in the process
        inline virtual void addEventContent() { InWarning( "Virtual method called" ); }
        /// Prepare the process for its integration over the whole phase space
        inline virtual void beforeComputeWeight() { Debugging( "Virtual method called" ); }
        /// Compute the weight for this point in the phase-space
        virtual double computeWeight() = 0;
        /// Fill the Event object with the particles' kinematics
        /// \param[in] symmetrise Symmetrise the event? (randomise the production of positively- and negatively-charged outgoing central particles)
        virtual void fillKinematics( bool symmetrise=false ) = 0;
        /// Return the number of dimensions on which the integration has to be performed
        /// \return Number of dimensions on which to integrate
        virtual unsigned int numDimensions( const Kinematics::ProcessMode& ) const = 0;
        /// Set the list of kinematic cuts to apply on the outgoing particles' final state
        /// \param[in] cuts_ The Cuts object containing the kinematic parameters
        inline virtual void setKinematics( const Kinematics& cuts ) { cuts_ = cuts; }

      public:
        /**
         * Sets the phase space point to compute the weight associated to it.
         * \brief Sets the phase space point to compute
         * \param[in] ndim The number of dimensions of the point in the phase space
         * \param[in] x[] The (@a ndim_)-dimensional point in the phase space on which the kinematics and the cross-section are computed
         */
        void setPoint( const unsigned int ndim, double* x );
        /// Dump the evaluated point's coordinates in the standard output stream
        void dumpPoint( const ExceptionType& et );
        /// Complete list of Particle with their role in the process for the point considered in the phase space, returned as an Event object.
        /// \return Event object containing all the generated Particle objects
        inline Event* event() { return event_; }

        ///Get the number of dimensions on which the integration is performed
        inline const unsigned int ndim() const { return num_dimensions_; }
        /// Get the value of a component of the @a fNumDimensions -dimensional point considered
        inline const double x( const unsigned int idx ) const {
          return ( idx>=num_dimensions_ ) ? -1. : x_[idx];
        }
        /// Get a human-readable name of the process considered
        inline const std::string& name() const { return name_; }
  
        /// Reset the total generation time and the number of events generated for this run
        inline void clearRun() {
          total_gen_time_ = 0.;
          num_gen_events_ = 0;
        }
        /// Add a new timing into the total generation time
        /// \param[in] gen_time Time to add (in seconds)
        inline void addGenerationTime( const float& gen_time ) {
          total_gen_time_ += gen_time;
          num_gen_events_++;
        }
        /// Return the total generation time for this run (in seconds)
        inline float totalGenerationTime() const { return total_gen_time_; }
        /// Total number of events already generated in this run
        inline unsigned int numGeneratedEvents() const { return num_gen_events_; }
  
      protected:
        /// Set the incoming and outgoing states to be defined in this process (and prepare the Event object accordingly)
        void setEventContent( const IncomingState& is, const OutgoingState& os );
        void formFactors( double q1, double q2, FormFactors& fp1, FormFactors& fp2 ) const;
 
        /// Get a list of pointers to the particles with a given role in the process
        /// \param[in] role role in the process for the particle to retrieve
        /// \return A vector of pointers to Particle objects associated to the role
        inline ParticlesRef particles( const Particle::Role& role ) { return event_->getByRole( role ); }
        /// Get the pointer to one particle in the event (using its role)
        inline Particle* particlePtr( const Particle::Role& role, unsigned int id=0 ) {
          if ( id==0 ) return event_->getOneByRole( role );
          ParticlesRef pp = event_->getByRole( role );
          if ( !pp.size() or id>pp.size() ) return 0;
          return pp.at( id );
        }
        /// Get the pointer to one particle in the event (using its identifier)
        inline Particle* particlePtr( unsigned int id ) { return event_->getById( id ); }

        // --- 
  
        /// Array of @a fNumDimensions components representing the point on which the weight in the cross-section is computed
        double* x_;
        /// List of incoming state particles (including intermediate partons)
        IncomingState incoming_state_;
        /// List of outgoing state particles
        OutgoingState outgoing_state_;
        /// \f$s\f$, squared centre of mass energy of the incoming particles' system, in \f$\mathrm{GeV}^2\f$
        double s_;
        /// \f$\sqrt s\f$, centre of mass energy of the incoming particles' system (in GeV)
        double sqs_;
        /// \f$m_1^2\f$, squared mass of the first proton-like incoming particle
        double w1_;
        /// \f$m_2^2\f$, squared mass of the second proton-like incoming particle
        double w2_;
        /// Virtuality of the first incoming photon
        double t1_;
        /// Virtuality of the second incoming photon
        double t2_;
        /// Invariant mass of the first proton-like outgoing particle (or remnant)
        double MX_;
        /// Invariant mass of the second proton-like outgoing particle (or remnant)
        double MY_;

        /// Number of dimensions on which the integration has to be performed.
        unsigned int num_dimensions_;
        /// Set of cuts to apply on the final phase space
        Kinematics cuts_;
        /// Event object containing all the information on the in- and outgoing particles
        Event* event_;
        /// Is the phase space point set?
        bool is_point_set_;
        /// Are the event's incoming particles set?
        bool is_incoming_state_set_;
        /// Are the event's outgoing particles set?
        bool is_outgoing_state_set_;
        /// Is the full event's kinematic set?
        bool is_kinematics_set_;
        /// Name of the process (useful for logging and debugging)
        std::string name_;
        /// Total generation time (in seconds)
        float total_gen_time_;
        /// Number of events already generated
        unsigned int num_gen_events_;

      private:
        /**
         * Is the system's kinematics well defined and compatible with the process ?
         * This check is mandatory to perform the (@a fNumDimensions)-dimensional point's cross-section computation.
         * \brief Is the system's kinematics well defined?
         * \return A boolean stating if the input kinematics and the final states are well-defined
         */
        inline bool isKinematicsDefined() {
          // check the incoming state
          if ( particles( Particle::IncomingBeam1 ).size()!=0 and particles( Particle::IncomingBeam2 ).size()!=0 ) {
            is_incoming_state_set_ = true;
          }
          // check the outgoing state
          if ( ( particles( Particle::OutgoingBeam1 ).size()!=0   and particles( Particle::OutgoingBeam2 ).size()!=0 )
           and ( particles( Particle::CentralParticle1 ).size()!=0 or particles( Particle::CentralParticle2 ).size()!=0) ) {
            is_outgoing_state_set_ = true;
          }
          // combine both states
          is_kinematics_set_ = is_incoming_state_set_ and is_outgoing_state_set_;
          return is_kinematics_set_;
        }
    };
  }
}

#endif
