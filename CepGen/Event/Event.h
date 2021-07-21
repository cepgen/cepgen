#ifndef CepGen_Event_Event_h
#define CepGen_Event_Event_h

#include "CepGen/Event/Particle.h"
#include <memory>

namespace cepgen {
  /**
   * Class containing all the information on the in- and outgoing particles' kinematics
   * \brief Kinematic information on the particles in the event
   */
  class Event {
  public:
    /// Build an empty event
    explicit Event(bool compressed = false);
    /// Copy constructor
    Event(const Event&) = default;
    /// Empty the whole event content
    void clear();
    /// Initialize an "empty" event collection
    void freeze();
    /// Restore the event to its "empty" state
    void restore();
    /// Is the event already without intermediate-channel information?
    bool compressed() const;
    /// Compress the event record
    Event compress() const;

    /// Callback function each event is required to pass through
    typedef std::function<void(const Event&, unsigned long long)> callback;

    /// Human-readable version of the event content
    friend std::ostream& operator<<(std::ostream&, const Event&);
    /// Dump all the known information on every Particle object contained in this Event container in the output stream
    void dump() const;
    /// Incoming beams centre-of-mass energy, in GeV
    double cmEnergy() const;

    //----- particles adders

    /// \brief Set the information on one particle in the process
    /// \param[in] part The Particle object to insert or modify in the event
    /// \param[in] replace Do we replace the particle if already present in the event or do we append another particle with the same role ?
    Particle& addParticle(Particle& part, bool replace = false);
    /// \brief Create a new particle in the event, with no kinematic information but the role it has to play in the process
    /// \param[in] role The role the particle will play in the process
    /// \param[in] replace Do we replace the particle if already present in the event or do we append another particle with the same role ?
    Particle& addParticle(Particle::Role role, bool replace = false);

    //----- particles retrievers

    /// Number of particles in the event
    size_t size() const;
    /// Vector of all particles in the event
    Particles particles() const;
    /// Vector of all stable particles in the event
    Particles stableParticles() const;
    /// List of references to Particle objects corresponding to a certain role in the process kinematics
    /// \param[in] role The role the particles have to play in the process
    Particles& operator[](Particle::Role role);
    /// Get a list of constant Particle objects corresponding to a certain role in the process kinematics
    const Particles& operator[](Particle::Role role) const;
    /// Get a list of particle identifiers in Event corresponding to a certain role in the process kinematics
    ParticlesIds ids(Particle::Role role) const;
    /// First Particle object with a given role in the event
    /// \param[in] role The role the particle has to play in the event
    Particle& oneWithRole(Particle::Role role);
    /// First constant Particle object with a given role in the event
    const Particle& oneWithRole(Particle::Role role) const;
    /// Reference to the Particle object corresponding to a unique identifier in the event
    /// \param[in] id The unique identifier to this particle in the event
    Particle& operator[](int id);
    /// Constant Particle reference object using its unique identifier
    /// \param[in] id Unique identifier of the particle in the event
    const Particle& operator[](int id) const;
    /// References to the Particle objects corresponding to the unique identifiers in the event
    /// \param[in] ids_ The unique identifiers to the particles to be selected in the event
    Particles operator[](const ParticlesIds& ids_) const;

    //----- general particles information retriever

    /// List of all parent Particle object for this given particle
    /// \param[in] part The particle for which the mother particles have to be retrieved
    Particles mothers(const Particle& part) const;
    /// List of all the daughters from a particle
    /// \param[in] part The particle for which the daughter particles have to be retrieved
    Particles daughters(const Particle& part) const;
    /// List of roles defined for the given event (really process-dependant for the central system)
    ParticleRoles roles() const;

    /// Number of trials before the event was "correctly" hadronised
    unsigned short num_hadronisation_trials;
    /// Time needed to generate the event at parton level (in seconds)
    float time_generation;
    /// Time needed to generate the hadronised (if needed) event (in seconds)
    float time_total;
    /// Event weight
    float weight;

  private:
    static constexpr double MIN_PRECISION = 1.e-10;
    /// Check if the event kinematics is properly defined
    void checkKinematics() const;
    /// List of particles in the event, mapped to their role in the process
    ParticlesMap particles_;
    /// Typical event indices structure
    struct NumParticles {
      size_t cs{0};   ///< Index of the first central system particle
      size_t op1{0};  ///< Index of the first positive-z outgoing beam state
      size_t op2{0};  ///< Index of the first negative-z outgoing beam state
    } evtcontent_{};
    /// Is the event "compressed"?
    bool compressed_;
  };
}  // namespace cepgen

#endif
