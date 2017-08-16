#ifndef CepGen_Physics_Event_h
#define CepGen_Physics_Event_h

#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <algorithm>

#include "Particle.h"

namespace CepGen
{
  /**
   * Class containing all the information on the in- and outgoing particles' kinematics
   * \brief Kinematic information on the particles in the event
   */
  class Event {
    public:
      Event();
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
      void init();
      /// Restore the event to its "empty" state
      void restore();
      /**
       * Returns the list of pointers to the Particle objects corresponding to a certain role in the process kinematics
       * \brief Gets a list of particles by their role in the event
       * \param[in] role_ The role the particles have to play in the process
       * \return A vector of pointers to the requested Particle objects
       */
      ParticlesRef getByRole( const Particle::Role& role_);
      /**
       * Returns the first Particle object in the particles list whose role corresponds to the given argument
       * \param[in] role_ The role the particle has to play in the event
       * \return A Particle object corresponding to the first particle found in this event
       */
      inline Particle* getOneByRole( const Particle::Role& role_ ) {
        ParticlesMap::iterator it = particles_.find( role_ );
        if ( it!=particles_.end() ) return &( it->second );
        return 0;
      };
      /**
       * Returns the pointer to the Particle object corresponding to a unique identifier in the event
       * \brief Gets one particle by its unique identifier in the event
       * \param[in] id_ The unique identifier to this particle in the event
       * \return A pointer to the requested Particle object
       */
      Particle* getById( int id_ );
      /// Get a const Particle object using its unique identifier
      /// \param[in] id_ Unique identifier of the particle in the event
      /// \return Constant object to be retrieved
      const Particle getConstById( int id_ ) const;
      /**
       * Returns the pointers to the Particle objects corresponding to the unique identifiers in the event
       * \brief Gets a vector of particles by their unique identifier in the event
       * \param[in] ids_ The unique identifiers to the particles to be selected in the event
       * \return A vector of pointers to the requested Particle objects
       */
      inline ParticlesRef getByIds( const ParticlesIds& ids_ ) {
        ParticlesRef out;
        for ( ParticlesIds::const_iterator id=ids_.begin(); id!=ids_.end(); id++ ) {
          out.push_back( getById( *id ) );
        }
        return out;
      }
      /**
       * Returns the pointer to the mother particle of any given Particle object in this event
       * \param[in] part_ The pointer to the Particle object from which we want to extract the mother particle
       * \return A pointer to the mother Particle object
       */
      inline ParticlesRef mothers( Particle* part_ ) {
        ParticlesRef out;
        const ParticlesIds moth = part_->mothersIds();
        for ( ParticlesIds::const_iterator m=moth.begin(); m!=moth.end(); m++ ) {
          out.push_back( getById( *m ) );
        }
        return out;
      }
      /// Return all const objects representing the mother particles of a given particle
      /// \param[in] part_ Particle object for which the mothers are retrieved
      /// \return Vector of Particle mother objects
      inline Particles constMothers( const Particle* part_ ) const {
        Particles out;
        const ParticlesIds moth = part_->mothersIds();
        for ( ParticlesIds::const_iterator m=moth.begin(); m!=moth.end(); m++ ) {
          out.push_back( getConstById( *m ) );
        }
        return out;
      }
      /// Get a vector containing all the daughters from a particle
      /// \param[in] part_ The particle for which the daughter particles have to be retrieved
      /// \return Vector of Particle objects containing all the daughters' kinematic information
      inline ParticlesRef daughters( const Particle* part_ ) { return getByIds( part_->daughters() ); };
      /// Get a list of roles for the given event (really process-dependant for the central system)
      /// \return Vector of integers corresponding to all the roles the particles can play in the event
      ParticleRoles roles() const;
      /// Set the information on one particle in the process
      /**
       * \param[in] part_ The Particle object to insert or modify in the event
       * \param[in] replace_ Do we replace the particle if already present in the event or do we append another particle with the same role ?
       * \return
       *  * 1 if a new Particle object has been inserted in the event
       *  * 0 if an existing Particle object has been modified
       *  * -1 if the requested role to edit is undefined or incorrect
       */
      int addParticle( Particle part_, bool replace_=false );
      /// Create a new particle in the event, with no kinematic information but the role it has to play in the process
      /**
       * \param[in] role_ The role the particle will play in the process
       * \param[in] replace_ Do we replace the particle if already present in the event or do we append another particle with the same role ?
       * \return
       *  * 1 if a new Particle object has been inserted in the event
       *  * 0 if an existing Particle object has been modified
       *  * -1 if the requested role to edit is undefined or incorrect
       */
      int addParticle( const Particle::Role& role_, bool replace_=false );
      /// Dump all the known information on every Particle object contained in this Event container in the output stream
      /// \param[in] stable_ Do we only show the stable particles in this event?
      void dump( bool stable_=false ) const;
      /// Get a vector of all particles in the event
      /// \return Vector containing all the pointers to the Particle objects contained in the event
      ParticlesRef particles();
      /// Get a vector of all particles in the event as const objects
      /// \return Vector containing all the const pointers to the Particle objects contained in the event
      Particles constParticles() const;
      /// Get the list of references to const particles contained in the event
      ConstParticlesRef constParticlesRef() const;
      /// Get a vector of all stable particles in the event
      /// \return Vector containing all the pointers to the stable Particle objects contained in the event
      ParticlesRef stableParticles();
      /// Number of particles in the event
      /// \return Integer number of particles in the event
      inline unsigned int numParticles() const { return particles_.size(); };
      /// Number of trials before the event was "correctly" hadronised
      int num_hadronisation_trials;
      /// Time needed to generate the event at parton level (in seconds)
      float time_generation;
      /// Time needed to generate the hadronised (if needed) event (in seconds)
      float time_total;

    private:
      /// List of particles in the event, mapped to their role in the process
      ParticlesMap particles_;
      /// Last particle in an "empty" event
      ParticlesMap::iterator last_particle_;
  };
}

#endif
