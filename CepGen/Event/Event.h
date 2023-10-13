/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2023  Laurent Forthomme
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CepGen_Event_Event_h
#define CepGen_Event_Event_h

#include <memory>

#include "CepGen/Event/Particle.h"

namespace cepgen {
  /// Container for the information on the in- and outgoing particles' kinematics
  class Event {
  public:
    explicit Event(bool compressed = false);  ///< Build an empty event
    Event(const Event&);                      ///< Copy constructor

    Event& operator=(const Event&);  ///< Assignment operator

    /// Build a trivial event with the minimal information
    /// \param[in] num_out_particles produced particles multiplicity (excluding outgoing beam remnants)
    static Event minimal(size_t num_out_particles = 1);

    void clear();             ///< Empty the whole event content
    void freeze();            ///< Initialize an "empty" event collection
    void restore();           ///< Restore the event to its "empty" state
    bool compressed() const;  ///< Is the event already without intermediate-channel information?
    Event compress() const;   ///< Compress the event record

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
    ParticleRef addParticle(Particle& part, bool replace = false);
    /// \brief Create a new particle in the event, with no kinematic information but the role it has to play in the process
    /// \param[in] role The role the particle will play in the process
    /// \param[in] replace Do we replace the particle if already present in the event or do we append another particle with the same role ?
    ParticleRef addParticle(Particle::Role role, bool replace = false);

    //----- particles retrievers

    size_t size() const;                        ///< Number of particles in the event
    Particles particles() const;                ///< Vector of all particles in the event
    Particles stableParticles() const;          ///< Vector of all stable particles in the event
    ParticlesMap& map() { return particles_; }  ///< Internal particles map retrieval operator

    /// List of references to Particle objects corresponding to a certain role in the process kinematics
    /// \param[in] role The role the particles have to play in the process
    ParticlesRefs operator[](Particle::Role role);
    /// Get a list of constant Particle objects corresponding to a certain role in the process kinematics
    const Particles& operator()(Particle::Role role) const;
    /// Get a list of particle identifiers in Event corresponding to a certain role in the process kinematics
    ParticlesIds ids(Particle::Role role) const;
    /// Check whether a particle role is represented in this event
    bool hasRole(Particle::Role role) const { return particles_.count(role) != 0; }
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

    /// Compute the missing momentum for central particles in this event
    Momentum missingMomentum() const;

    //----- general particles information retriever

    /// List of all parent Particle object for this given particle
    /// \param[in] part The particle for which the mother particles have to be retrieved
    Particles mothers(const Particle& part) const;
    /// List of all the daughters from a particle
    /// \param[in] part The particle for which the daughter particles have to be retrieved
    Particles daughters(const Particle& part) const;
    /// List of roles defined for the given event (really process-dependant for the central system)
    ParticleRoles roles() const;

    unsigned short num_hadronisation_trials{0};  ///< Number of trials before the event was "correctly" hadronised
    float time_generation{-1.};                  ///< Time (in s) to generate the event at parton level
    float time_total{-1.};                       ///< Time (in s) to generate the (possibly modified/hadronised) event
    float weight{1.};                            ///< Event weight
    float alpha_em{constants::ALPHA_EM};         ///< Electromagnetic coupling constant
    float alpha_s{constants::ALPHA_QCD};         ///< Strong coupling constant

  private:
    static constexpr double MIN_PRECISION = 1.e-10;
    void checkKinematics() const;  ///< Check if the event kinematics is properly defined
    ParticlesMap particles_;       ///< List of particles in the event, mapped to their role in the process
    /// Typical event indices structure
    struct NumParticles {
      size_t cs{0};   ///< Index of the first central system particle
      size_t op1{0};  ///< Index of the first positive-z outgoing beam state
      size_t op2{0};  ///< Index of the first negative-z outgoing beam state
    } evtcontent_{};
    bool compressed_{false};  ///< Is the event "compressed"?
  };
}  // namespace cepgen

#endif
