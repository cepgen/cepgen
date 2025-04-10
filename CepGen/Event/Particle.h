/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2025  Laurent Forthomme
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

#ifndef CepGen_Event_Particle_h
#define CepGen_Event_Particle_h

#include <set>

#include "CepGen/Physics/Momentum.h"
#include "CepGen/Physics/ParticleProperties.h"
#include "CepGen/Utils/Hasher.h"

namespace cepgen {
  typedef std::set<int> ParticlesIds;  ///< A set of integer-type particle identifiers

  /// Kinematic information for one particle
  class Particle {
  public:
    /// Internal status code for a particle
    enum class Status {
      PrimordialIncoming = -9,  ///< Incoming beam particle
      DebugResonance = -5,      ///< Intermediate resonance (for processes developers)
      Resonance = -4,           ///< Already decayed intermediate resonance
      Fragmented = -3,          ///< Already fragmented outgoing beam
      Propagator = -2,          ///< Generic propagator
      Incoming = -1,            ///< Incoming parton
      Undefined = 0,            ///< Undefined particle
      FinalState = 1,           ///< Stable, final state particle
      Undecayed = 2,            ///< Particle to be decayed externally
      Unfragmented = 3          ///< Particle to be hadronised externally
    };
    friend std::ostream& operator<<(std::ostream& os, const Status&);  ///< Human-readable particle's status
    /// Role of the particle in the process
    enum class Role {
      UnknownRole = -1,   ///< Undefined role
      IncomingBeam1 = 1,  ///< \f$z>0\f$ incoming beam particle
      IncomingBeam2 = 2,  ///< \f$z<0\f$ incoming beam particle
      OutgoingBeam1 = 3,  ///< \f$z<0\f$ outgoing beam state/particle
      OutgoingBeam2 = 5,  ///< \f$z>0\f$ outgoing beam state/particle
      CentralSystem = 6,  ///< Central particles system
      Intermediate = 4,   ///< Intermediate two-parton system
      Parton1 = 41,       ///< \f$z>0\f$ beam incoming parton
      Parton2 = 42        ///< \f$z<0\f$ beam incoming parton
    };
    friend std::ostream& operator<<(std::ostream& os, const Role&);  ///< Human-readable particle's role in the event

    //----- static getters

    /// Build using the role of the particle in the process and its PDG id
    /// \param[in] role Role of the particle in the process
    /// \param[in] id PDG identifier
    /// \param[in] st Current status
    explicit Particle(Role role = Role::UnknownRole, pdgid_t id = 0, Status st = Status::Undefined);

    bool operator<(const Particle&) const;   ///< Comparison operator (from unique identifier)
    bool operator==(const Particle&) const;  ///< Equality operator

    // --- general particle properties

    inline int id() const { return id_; }  ///< Unique identifier (in a Event object context)
    /// Set the particle unique identifier in an event
    inline Particle& setId(int id) {
      id_ = id;
      return *this;
    }
    float charge() const;  ///< Electric charge (given as a float number, for the quarks and bound states)
    /// Set whether we are coping with the particle or its antiparticle
    inline Particle& setAntiparticle(bool anti) {
      antiparticle_ = anti;
      return *this;
    }
    inline Role role() const { return role_; }  ///< Role in the considered process
    /// Set the particle role in the process
    inline Particle& setRole(const Role& role) {
      role_ = role;
      return *this;
    }
    inline Status status() const { return static_cast<Status>(status_); }  ///< Particle status
    /// Set the particle decay/stability status
    inline Particle& setStatus(Status status) {
      status_ = static_cast<int>(status);
      return *this;
    }
    /// Set the particle decay/stability status
    inline Particle& setStatus(int status) {
      status_ = status;
      return *this;
    }

    /// Set the PDG identifier (along with the particle's electric charge)
    /// \param[in] pdg PDG identifier
    /// \param[in] ch Electric charge (0, 1, or -1)
    Particle& setPdgId(pdgid_t pdg, short ch = 0);
    pdgid_t pdgId() const;  ///< Retrieve the objectified PDG identifier
    /// Set the PDG identifier (along with the particle's electric charge)
    /// \param[in] pdg_id PDG identifier (incl. electric charge in e)
    Particle& setIntegerPdgId(long pdg_id);
    long integerPdgId() const;  ///< Retrieve the integer value of the PDG identifier

    inline float helicity() const { return helicity_; }  ///< Particle's helicity
    /// Set the helicity of the particle
    inline Particle& setHelicity(float heli) {
      helicity_ = heli;
      return *this;
    }

    inline Momentum& momentum() { return momentum_; }  ///< Retrieve the momentum object associated with this particle
    /// Retrieve the momentum object associated with this particle
    inline const Momentum& momentum() const { return momentum_; }
    /// Associate a momentum object to this particle
    /// \param[in] off_shell allow the 4-momentum mass to compensate for E-p balance?
    Particle& setMomentum(const Momentum&, bool off_shell = false);
    /// Set the 3- or 4-momentum associated to the particle
    /// \param[in] px Momentum along the \f$x\f$-axis, in GeV/c
    /// \param[in] py Momentum along the \f$y\f$-axis, in GeV/c
    /// \param[in] pz Momentum along the \f$z\f$-axis, in GeV/c
    /// \param[in] energy Energy, in GeV
    Particle& setMomentum(double px, double py, double pz, double energy = -1.);
    /// Set the 4-momentum associated to the particle
    /// \param[in] p 4-momentum
    inline Particle& setMomentum(double p[4]) { return setMomentum(p[0], p[1], p[2], p[3]); }
    bool valid() const;  ///< Is this particle a valid particle which can be used for kinematic computations?

    // --- particle relations

    inline bool primary() const { return mothers_.empty(); }    ///< Is this particle a primary particle?
    Particle& addMother(Particle& mother_particle);             ///< Set the mother particle
    inline ParticlesIds mothers() const { return mothers_; }    ///< Identifier to the mother particles
    inline ParticlesIds& mothers() { return mothers_; }         ///< Identifier to the mother particles
    Particle& addChild(Particle& child_particle);               ///< Add a decay product
    inline ParticlesIds children() const { return children_; }  ///< Identifiers list of all child particles
    inline ParticlesIds& children() { return children_; }       ///< Identifiers list of all child particles

    // --- global particle information extraction

    friend std::ostream& operator<<(std::ostream&, const Particle&);  ///< Human-readable dump of particle information

  protected:
    int id_{-1};                                       ///< Unique identifier in an event
    bool antiparticle_{false};                         ///< Are we dealing with the particle or antiparticle?
    Momentum momentum_;                                ///< Momentum properties handler
    float helicity_{0.};                               ///< Helicity
    Role role_{Role::UnknownRole};                     ///< Role in the process
    int status_{static_cast<int>(Status::Undefined)};  ///< Decay/stability status
    ParticlesIds mothers_{};                           ///< List of mother particles
    ParticlesIds children_{};                          ///< List of child particles
    pdgid_t pdg_id_{0};                                ///< PDG id
  };

  // --- particle containers

  typedef std::reference_wrapper<Particle> ParticleRef;  ///< Reference to a Particle object
  typedef std::vector<Particle> Particles;               ///< List of Particle objects
  typedef std::vector<ParticleRef> ParticlesRefs;        ///< List of references to Particle objects
  typedef std::vector<Particle::Role> ParticleRoles;     ///< List of particles' roles

  /// Map between a particle's role and its associated Particle object
  class ParticlesMap : public std::unordered_map<Particle::Role, Particles, utils::EnumHash<Particle::Role> > {
  public:
    ParticlesMap() = default;
    ParticlesMap(const ParticlesMap&);             ///< Copy constructor
    ParticlesMap& operator=(const ParticlesMap&);  ///< Assignment operator
    ~ParticlesMap() = default;
  };
}  // namespace cepgen

#endif
