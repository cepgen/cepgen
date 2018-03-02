#ifndef CepGen_Event_Event_h
#define CepGen_Event_Event_h

#include "Particle.h"
#include "CepGen/Core/Logger.h"

namespace CepGen
{
  /**
   * Class containing all the information on the in- and outgoing particles' kinematics
   * \brief Kinematic information on the particles in the event
   */
  class Event {
    public:
      Event();
      Event( const Event& );
      ~Event();
      /**
       * \brief Copies all the relevant quantities from one Event object to another
       */
      Event& operator=( const Event& );
      /**
       * \brief Empties the whole event content
       */
      void clear();
      /// Initialize an "empty" event collection
      void freeze();
      /// Restore the event to its "empty" state
      void restore();

      /// Dump all the known information on every Particle object contained in this Event container in the output stream
      /// \param[out] os Output stream where to dump the information
      /// \param[in] stable_ Do we only show the stable particles in this event?
      void dump( std::ostream& os = Logger::get().outputStream, bool stable_=false ) const;

      double cmEnergy() const { return CMEnergy( getOneByRole( Particle::IncomingBeam1 ), getOneByRole( Particle::IncomingBeam2 ) ); }

      //----- particles adders

      /// Set the information on one particle in the process
      /// \param[in] part The Particle object to insert or modify in the event
      /// \param[in] replace Do we replace the particle if already present in the event or do we append another particle with the same role ?
      Particle& addParticle( Particle& part, bool replace=false );
      /// \brief Create a new particle in the event, with no kinematic information but the role it has to play in the process
      /// \param[in] role The role the particle will play in the process
      /// \param[in] replace Do we replace the particle if already present in the event or do we append another particle with the same role ?
      Particle& addParticle( Particle::Role role, bool replace=false );

      //----- particles retrievers

      /// Number of particles in the event
      size_t numParticles() const;
      /// \brief Vector of all particles in the event
      const Particles particles() const;
      /// \brief Vector of all stable particles in the event
      const Particles stableParticles() const;
      /**
       * Returns the list of Particle objects corresponding to a certain role in the process kinematics
       * \brief Gets a list of particles by their role in the event
       * \param[in] role The role the particles have to play in the process
       * \return A vector of references to the requested Particle objects
       */
      Particles& getByRole( Particle::Role role );
      const Particles& getByRole( Particle::Role role ) const;
      ParticlesIds getIdsByRole( Particle::Role role ) const;
      /**
       * Returns the first Particle object in the particles list whose role corresponds to the given argument
       * \param[in] role The role the particle has to play in the event
       * \return A Particle object corresponding to the first particle with the role
       */
      Particle& getOneByRole( Particle::Role role );
      const Particle& getOneByRole( Particle::Role role ) const;
      /**
       * Returns the reference to the Particle object corresponding to a unique identifier in the event
       * \brief Gets one particle by its unique identifier in the event
       * \param[in] id_ The unique identifier to this particle in the event
       * \return A reference to the requested Particle object
       */
      Particle& getById( int id_ );
      /// Get a const Particle object using its unique identifier
      /// \param[in] id_ Unique identifier of the particle in the event
      /// \return Constant object to be retrieved
      const Particle& getConstById( int id_ ) const;
      /**
       * Returns the references to the Particle objects corresponding to the unique identifiers in the event
       * \brief Gets a vector of particles by their unique identifier in the event
       * \param[in] ids_ The unique identifiers to the particles to be selected in the event
       * \return A vector of references to the requested Particle objects
       */
      Particles getByIds( const ParticlesIds& ids_ ) const;
      /**
       * Returns the list of mother particles of any given Particle object in this event
       * \param[in] part The reference to the Particle object from which we want to extract the mother particles
       * \return A list of parenting Particle object
       */

      //----- general particles information retriever

      Particles mothers( const Particle& part );
      /// Get a vector containing all the daughters from a particle
      /// \param[in] part The particle for which the daughter particles have to be retrieved
      /// \return Vector of Particle objects containing all the daughters' kinematic information
      Particles daughters( const Particle& part );
      /// Get a list of roles for the given event (really process-dependant for the central system)
      /// \return Vector of integers corresponding to all the roles the particles can play in the event
      ParticleRoles roles() const;

      /// Number of trials before the event was "correctly" hadronised
      unsigned short num_hadronisation_trials;
      /// Time needed to generate the event at parton level (in seconds)
      float time_generation;
      /// Time needed to generate the hadronised (if needed) event (in seconds)
      float time_total;

    private:
      static constexpr double minimal_precision_ = 1.e-10;
      /// Check if the event kinematics is properly defined
      void checkKinematics() const;
      /// List of particles in the event, mapped to their role in the process
      ParticlesMap particles_;
      /// Last particle in an "empty" event
      struct num_particles {
        num_particles() : cs( 0 ), op1( 0 ), op2( 0 ) {}
        num_particles( const num_particles& np ) : cs( np.cs ), op1( np.op1 ), op2( np.op2 ) {}
        unsigned short cs, op1, op2;
      };
      num_particles evtcontent_;
  };
}

#endif
