#ifndef Event_h
#define Event_h

#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <algorithm>

#include "Particle.h"

/**
 * Class containing all the information on the in- and outgoing particles' kinematics
 * @brief Kinematic information on the particles in the event
 */
class Event {
  public:
    Event();
    ~Event();
    /**
     * @brief Copies all the relevant quantities from one Event object to another
     */
    Event& operator=(const Event&);
    /**
     * @brief Empties the whole event content
     */
    void clear();
    /// Initialize an "empty" event collection
    void Init();
    /// Restore the event to its "empty" state
    void Restore();
    /**
     * Returns the list of pointers to the Particle objects corresponding to a certain role in the process kinematics
     * @brief Gets a list of particles by their role in the event
     * @param[in] role_ The role the particles have to play in the process
     * @return A vector of pointers to the requested Particle objects
     */
    ParticlesRef GetByRole(Particle::Role role_);
    /**
     * Returns the first Particle object in the particles list whose role corresponds to the given argument
     * @param[in] role_ The role the particle has to play in the event
     * @return A Particle object corresponding to the first particle found in this event
     */
    inline Particle* GetOneByRole(Particle::Role role_) {
      ParticlesMap::iterator it = fParticles.find(role_);
      if (it!=fParticles.end()) return &(it->second);
      return 0;
    };
    inline const Particle* GetOneByRole(Particle::Role role_) const {
      ParticlesMap::const_iterator it = fParticles.find(role_);
      if (it!=fParticles.end()) return const_cast<const Particle*>(&(it->second));
      return 0;
    };
    /**
     * Returns the pointer to the Particle object corresponding to a unique identifier in the event
     * @brief Gets one particle by its unique identifier in the event
     * @param[in] id_ The unique identifier to this particle in the event
     * @return A pointer to the requested Particle object
     */
    Particle* GetById(int id_);
    /// Get a const Particle object using its unique identifier
    /// \param[in] id_ Unique identifier of the particle in the event
    /// \return Constant object to be retrieved
    const Particle GetConstById(int id_) const;
    /**
     * Returns the pointers to the Particle objects corresponding to the unique identifiers in the event
     * @brief Gets a vector of particles by their unique identifier in the event
     * @param[in] ids_ The unique identifiers to the particles to be selected in the event
     * @return A vector of pointers to the requested Particle objects
     */
    inline ParticlesRef GetByIds(std::vector<int> ids_) {
      std::vector<int>::iterator id;
      ParticlesRef out;
      for (id=ids_.begin(); id!=ids_.end(); id++) out.push_back(GetById(*id));
      return out;
    }
    /**
     * Returns the pointer to the mother particle of any given Particle object in this event
     * @param[in] part_ The pointer to the Particle object from which we want to extract the mother particle
     * @return A pointer to the mother Particle object
     */
    inline ParticlesRef GetMothers(Particle* part_) {
      ParticlesRef out;
      const ParticlesIds moth = part_->GetMothersIds();
      ParticlesIds::iterator m;
      for (m=moth.begin(); m!=moth.end(); m++) {
      	out.push_back(GetById(*m));
      }
      return out;
    }; // FIXME
    /// Return all const objects representing the mother particles of a given particle
    /// \param[in] part_ Particle object for which the mothers are retrieved
    /// \return Vector of Particle mother objects
    inline Particles GetConstMothers(Particle* part_) const {
      Particles out;
      const ParticlesIds moth = part_->GetMothersIds();
      ParticlesIds::iterator m;
      for (m=moth.begin(); m!=moth.end(); m++) {
      	out.push_back(GetConstById(*m));
      }
      return out;
    }; // FIXME
    /// Get a vector containing all the daughters from a particle
    /// \param[in] part_ The particle for which the daughter particles have to be retrieved
    /// \return Vector of Particle objects containing all the daughters' kinematic information
    inline ParticlesRef GetDaughters(Particle* part_) { return GetByIds(part_->GetDaughters()); };
    /// Get a list of roles for the given event (really process-dependant for the central system)
    /// \return Vector of integers corresponding to all the roles the particles can play in the event
    ParticleRoles GetRoles() const;
    /// Set the information on one particle in the process
    /**
     * \param[in] part_ The Particle object to insert or modify in the event
     * \param[in] replace_ Do we replace the particle if already present in the event or do we append another particle with the same role ?
     * \return
     *  * 1 if a new Particle object has been inserted in the event
     *  * 0 if an existing Particle object has been modified
     *  * -1 if the requested role to edit is undefined or incorrect
     */
    int AddParticle(Particle part_, bool replace_=false);
    /// Create a new particle in the event, with no kinematic information but the role it has to play in the process
    /**
     * \param[in] role_ The role the particle will play in the process
     * \param[in] replace_ Do we replace the particle if already present in the event or do we append another particle with the same role ?
     * \return
     *  * 1 if a new Particle object has been inserted in the event
     *  * 0 if an existing Particle object has been modified
     *  * -1 if the requested role to edit is undefined or incorrect
     */
    int AddParticle(Particle::Role role_, bool replace_=false);
    //HEPEUP GetHEPEUP() const;
    /// Store in a file (raw format) all the kinematics on the outgoing leptons
    /// \param[in] weight_ Weight of the event
    void Store(std::ofstream*,double weight_=1.);
    //void Hadronise(std::string algo_="");
    /// Dump all the known information on every Particle object contained in this Event container in the output stream
    /// \param[in] stable_ Do we only show the stable particles in this event?
    void Dump(bool stable_=false) const;
    /// Get a vector of all particles in the event
    /// \return Vector containing all the pointers to the Particle objects contained in the event
    ParticlesRef GetParticles();
    /// Get a vector of all particles in the event as const objects
    /// \return Vector containing all the const pointers to the Particle objects contained in the event
    Particles GetConstParticles() const;
    ConstParticlesRef GetConstParticlesRef() const;
    /// Get a vector of all stable particles in the event
    /// \return Vector containing all the pointers to the stable Particle objects contained in the event
    ParticlesRef GetStableParticles();
    /// Number of particles in the event
    /// \return Integer number of particles in the event
    inline unsigned int NumParticles() const { return fParticles.size(); };
    /// Number of trials before the event was "correctly" hadronised
    int num_hadronisation_trials;
    /// Time needed to generate the event at parton level (in seconds)
    float time_generation;
    /// Time needed to generate the hadronised (if needed) event (in seconds)
    float time_total;
    //HEPEUP event_info;
  private:
    /// List of particles in the event, mapped to their role in the process
    ParticlesMap fParticles;
    /// Last particle in an "empty" event
    ParticlesMap::iterator fLastParticle;
    /// Empty particle returned to the get-ers if no particle matches the requirements
    Particle* np;
};

#endif
