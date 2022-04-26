/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
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
#include <unordered_map>
#include <vector>

#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/Momentum.h"
#include "CepGen/Physics/ParticleProperties.h"
#include "CepGen/Utils/Hasher.h"

namespace cepgen {
  /// A set of integer-type particle identifiers
  typedef std::set<int> ParticlesIds;

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
    /// Role of the particle in the process
    enum Role {
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
    /// Human-readable format for a particle's role in the event
    friend std::ostream& operator<<(std::ostream& os, const Role& rl);

    //----- static getters

    /// Convert a polar angle to a pseudo-rapidity
    static double thetaToEta(double theta);
    /// Convert a pseudo-rapidity to a polar angle
    static double etaToTheta(double eta);
    /// Convert a pseudo-rapidity to a rapidity
    static double etaToY(double eta_, double m_, double pt_);

    Particle() = default;
    /// Build using the role of the particle in the process and its PDG id
    /// \param[in] role Role of the particle in the process
    /// \param[in] id PDG identifier
    /// \param[in] st Current status
    Particle(Role role, pdgid_t id, Status st = Status::Undefined);
    /// Copy constructor
    Particle(const Particle&);
    inline ~Particle() = default;
    Particle& operator=(const Particle&) = default;  ///< Assignment operator
    /// Comparison operator (from unique identifier)
    bool operator<(const Particle& rhs) const;
    /// Comparison operator (from their reference's unique identifier)
    //bool operator<( Particle *rhs ) const { return ( id < rhs->id ); }

    // --- general particle properties

    /// Unique identifier (in a Event object context)
    int id() const { return id_; }
    //void setId( int id ) { id_ = id; }
    /// Set the particle unique identifier in an event
    Particle& setId(int id) {
      id_ = id;
      return *this;
    }
    /// Electric charge (given as a float number, for the quarks and bound states)
    float charge() const;
    /// Set the electric charge sign (+-1 for charged or 0 for neutral particles)
    Particle& setChargeSign(int sign) {
      charge_sign_ = sign;
      return *this;
    }
    /// Role in the considered process
    Role role() const { return role_; }
    /// Set the particle role in the process
    Particle& setRole(const Role& role) {
      role_ = role;
      return *this;
    }
    /**
       * Codes 1-10 correspond to currently existing partons/particles, and larger codes contain partons/particles which no longer exist, or other kinds of event information
       * \brief Particle status
       */
    Status status() const { return (Status)status_; }
    /// Set the particle decay/stability status
    Particle& setStatus(Status status) {
      status_ = (int)status;
      return *this;
    }
    /// Set the particle decay/stability status
    Particle& setStatus(int status) {
      status_ = status;
      return *this;
    }

    /// Set the PDG identifier (along with the particle's electric charge)
    /// \param[in] pdg PDG identifier
    /// \param[in] ch Electric charge (0, 1, or -1)
    Particle& setPdgId(pdgid_t pdg, short ch = 0);
    /// Set the PDG identifier (along with the particle's electric charge)
    /// \param[in] pdg_id PDG identifier (incl. electric charge in e)
    Particle& setPdgId(long pdg_id);
    /// Retrieve the objectified PDG identifier
    pdgid_t pdgId() const;
    /// Retrieve the integer value of the PDG identifier
    int integerPdgId() const;
    /// Particle's helicity
    float helicity() const { return helicity_; }
    /// Set the helicity of the particle
    Particle& setHelicity(float heli) {
      helicity_ = heli;
      return *this;
    }
    /// Particle mass in GeV/c\f$^2\f$
    /// \return Particle's mass
    inline double mass() const { return mass_; };
    /// Compute the particle mass
    /// \param[in] off_shell Allow the particle to be produced off-shell?
    /// \note This method ensures that the kinematics is properly set (the mass is set according to the energy and the momentum in priority)
    Particle& computeMass(bool off_shell = false);
    /// Set the particle mass, in GeV/c\f$^2\f$
    /// \param m Mass in GeV/c\f$^2\f$
    /// \note This method ensures that the kinematics is properly set (the mass is set according to the energy and the momentum in priority)
    Particle& setMass(double m = -1.);
    /// Particle squared mass, in GeV\f$^2\f$/c\f$^4\f$
    inline double mass2() const { return mass_ * mass_; };
    /// Retrieve the momentum object associated with this particle
    inline Momentum& momentum() { return momentum_; }
    /// Retrieve the momentum object associated with this particle
    inline Momentum momentum() const { return momentum_; }
    /// Associate a momentum object to this particle
    Particle& setMomentum(const Momentum& mom, bool offshell = false);
    /**
       * \brief Set the 3- or 4-momentum associated to the particle
       * \param[in] px Momentum along the \f$x\f$-axis, in GeV/c
       * \param[in] py Momentum along the \f$y\f$-axis, in GeV/c
       * \param[in] pz Momentum along the \f$z\f$-axis, in GeV/c
       * \param[in] e Energy, in GeV
       */
    Particle& setMomentum(double px, double py, double pz, double e = -1.);
    /// Set the 4-momentum associated to the particle
    /// \param[in] p 4-momentum
    inline Particle& setMomentum(double p[4]) { return setMomentum(p[0], p[1], p[2], p[3]); }
    /// Set the particle's energy
    /// \param[in] e Energy, in GeV
    Particle& setEnergy(double e = -1.);
    /// Get the particle's energy, in GeV
    double energy() const;
    /// Get the particle's squared energy, in GeV\f$^2\f$
    inline double energy2() const { return energy() * energy(); };
    /// Is this particle a valid particle which can be used for kinematic computations?
    bool valid();

    // --- particle relations

    /// Is this particle a primary particle?
    inline bool primary() const { return mothers_.empty(); }
    /// Clear the particle parentage
    Particle& clearMothers();
    /// Set the mother particle
    /// \param[in] part A Particle object containing all the information on the mother particle
    Particle& addMother(Particle& part);
    /// Get the unique identifier to the mother particle from which this particle arises
    /// \return An integer representing the unique identifier to the mother of this particle in the event
    inline ParticlesIds mothers() const { return mothers_; }
    /// Remove the decay products linking
    Particle& clearDaughters();
    /**
       * \brief Add a decay product
       * \param[in] part The Particle object in which this particle will disintegrate or convert
       * \return A boolean stating if the particle has been added to the daughters list or if it was already present before
       */
    Particle& addDaughter(Particle& part);
    /// Gets the number of daughter particles
    inline size_t numDaughters() const { return daughters_.size(); };
    /// Get an identifiers list all daughter particles
    /// \return An integer vector containing all the daughters' unique identifier in the event
    inline ParticlesIds daughters() const { return daughters_; }

    // --- global particle information extraction

    /// Dump all the information on this particle into the standard output stream
    friend std::ostream& operator<<(std::ostream&, const Particle&);

  protected:
    /// Unique identifier in an event
    int id_{-1};
    /// Electric charge (+-1 or 0)
    short charge_sign_{1};
    /// Momentum properties handler
    Momentum momentum_;
    /// Mass, in GeV/c\f$^2\f$
    double mass_{-1.};
    /// Helicity
    float helicity_{0.};
    /// Role in the process
    Role role_{UnknownRole};
    /// Decay/stability status
    int status_{(int)Status::Undefined};
    /// List of mother particles
    ParticlesIds mothers_;
    /// List of daughter particles
    ParticlesIds daughters_;
    /// PDG id
    pdgid_t pdg_id_{(pdgid_t)0};
    /// Collection of standard, bare-level physical properties
    ParticleProperties phys_prop_;
  };

  /// Compute the centre of mass energy of two particles (incoming or outgoing states)
  double CMEnergy(const Particle& p1, const Particle& p2);

  //bool operator<( const Particle& a, const Particle& b ) { return a.id<b.id; }

  // --- particle containers

  /// List of Particle objects
  typedef std::vector<Particle> Particles;
  /// List of references to Particle objects
  typedef std::vector<std::reference_wrapper<Particle> > ParticlesRefs;
  /// List of particles' roles
  typedef std::vector<Particle::Role> ParticleRoles;
  /// Map between a particle's role and its associated Particle object
  typedef std::unordered_map<Particle::Role, Particles, utils::EnumHash<Particle::Role> > ParticlesMap;
}  // namespace cepgen

#endif
