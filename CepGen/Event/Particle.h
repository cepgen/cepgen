#ifndef CepGen_Event_Particle_h
#define CepGen_Event_Particle_h

#include "CepGen/Physics/ParticleProperties.h"

#include <set>
#include <unordered_map>
#include <vector>

namespace cepgen
{

  /// A set of integer-type particle identifiers
  typedef std::set<int> ParticlesIds;

  /// Kinematic information for one particle
  class Particle {
    public:
      /// Internal status code for a particle
      enum class Status {
        PrimordialIncoming = -9, ///< Incoming beam particle
        DebugResonance = -5, ///< Intermediate resonance (for processes developers)
        Resonance = -4, ///< Already decayed intermediate resonance
        Fragmented = -3, ///< Already fragmented outgoing beam
        Propagator = -2, ///< Generic propagator
        Incoming = -1, ///< Incoming parton
        Undefined = 0, ///< Undefined particle
        FinalState = 1, ///< Stable, final state particle
        Undecayed = 2, ///< Particle to be decayed externally
        Unfragmented = 3 ///< Particle to be hadronised externally
      };
      /// Role of the particle in the process
      enum Role {
        UnknownRole = -1, ///< Undefined role
        IncomingBeam1 = 1, ///< \f$z>0\f$ incoming beam particle
        IncomingBeam2 = 2, ///< \f$z<0\f$ incoming beam particle
        OutgoingBeam1 = 3, ///< \f$z<0\f$ outgoing beam state/particle
        OutgoingBeam2 = 5, ///< \f$z>0\f$ outgoing beam state/particle
        CentralSystem = 6, ///< Central particles system
        Intermediate = 4, ///< Intermediate two-parton system
        Parton1 = 41, ///< \f$z>0\f$ beam incoming parton
        Parton2 = 42 ///< \f$z<0\f$ beam incoming parton
      };
      /**
       * Container for a particle's 4-momentum, along with useful methods to ease the development of any matrix element level generator
       * \brief 4-momentum for a particle
       * \date Dec 2015
       * \author Laurent Forthomme <laurent.forthomme@cern.ch>
       */
      class Momentum {
        public:
          /// Build a 4-momentum at rest with an invalid energy (no mass information known)
          Momentum();
          /// Build a 4-momentum using its 3-momentum coordinates and its energy
          Momentum( double x, double y, double z, double t = -1. );
          /// Build a 4-momentum using its 3-momentum coordinates and its energy
          Momentum( double* p );

          // --- static definitions

          /// Build a 3-momentum from its three pseudo-cylindric coordinates
          static Momentum fromPtEtaPhi( double pt, double eta, double phi, double e = -1. );
          /// Build a 4-momentum from its scalar momentum, and its polar and azimuthal angles
          static Momentum fromPThetaPhi( double p, double theta, double phi, double e = -1. );
          /// Build a 4-momentum from its four momentum and energy coordinates
          static Momentum fromPxPyPzE( double px, double py, double pz, double e );
          /// Build a 4-momentum from its three momentum coordinates and mass
          static Momentum fromPxPyPzM( double px, double py, double pz, double m );
          /// Build a 4-momentum from its transverse momentum, rapidity and mass
          static Momentum fromPxPyYM( double px, double py, double rap, double m );

          // --- vector and scalar operators

          /// Scalar product of the 3-momentum with another 3-momentum
          double threeProduct( const Momentum& ) const;
          /// Scalar product of the 4-momentum with another 4-momentum
          double fourProduct( const Momentum& ) const;
          /// Vector product of the 3-momentum with another 3-momentum
          double crossProduct( const Momentum& ) const;
          /// Add a 4-momentum through a 4-vector sum
          Momentum& operator+=( const Momentum& );
          /// Subtract a 4-momentum through a 4-vector sum
          Momentum& operator-=( const Momentum& );
          /// Scalar product of the 3-momentum with another 3-momentum
          double operator*=( const Momentum& );
          /// Multiply all 4-momentum coordinates by a scalar
          Momentum& operator*=( double c );
          /// Equality operator
          bool operator==( const Momentum& ) const;
          /// Human-readable format for a particle's momentum
          friend std::ostream& operator<<( std::ostream& os, const Particle::Momentum& mom );

          Momentum& betaGammaBoost( double gamma, double betagamma );
          /// Forward Lorentz boost
          Momentum& lorentzBoost( const Particle::Momentum& p );

          // --- setters and getters

          /// Set all the components of the 4-momentum (in GeV)
          void setP( double px, double py, double pz, double e );
          /// Set all the components of the 3-momentum (in GeV)
          void setP( double px, double py, double pz );
          /// Set the energy (in GeV)
          inline void setEnergy( double e ) { energy_ = e; }
          /// Compute the energy from the mass
          inline void setMass( double m ) { setMass2( m*m ); }
          /// Compute the energy from the mass
          void setMass2( double m2 );
          /// Get one component of the 4-momentum (in GeV)
          double operator[]( const unsigned int i ) const;
          /// Get one component of the 4-momentum (in GeV)
          double& operator[]( const unsigned int i );
          /// Momentum along the \f$x\f$-axis (in GeV)
          inline double px() const { return px_; }
          /// Momentum along the \f$y\f$-axis (in GeV)
          inline double py() const { return py_; }
          /// Longitudinal momentum (in GeV)
          inline double pz() const { return pz_; }
          /// Transverse momentum (in GeV)
          double pt() const;
          /// Squared transverse momentum (in GeV\f$^2\f$)
          double pt2() const;
          /// 4-vector of double precision floats (in GeV)
          const std::vector<double> pVector() const;
          /// 3-momentum norm (in GeV)
          inline double p() const { return p_; }
          /// Squared 3-momentum norm (in GeV\f$^2\f$)
          inline double p2() const { return p_*p_; }
          /// Energy (in GeV)
          inline double energy() const { return energy_; }
          /// Squared energy (in GeV\f$^2\f$)
          inline double energy2() const { return energy_*energy_; }
          /// Squared mass (in GeV\f$^2\f$) as computed from its energy and momentum
          inline double mass2() const { return energy2()-p2(); }
          /// Mass (in GeV) as computed from its energy and momentum
          /// \note Returns \f$-\sqrt{|E^2-\mathbf{p}^2|}<0\f$ if \f$\mathbf{p}^2>E^2\f$
          double mass() const;
          /// Polar angle (angle with respect to the longitudinal direction)
          double theta() const;
          /// Azimutal angle (angle in the transverse plane)
          double phi() const;
          /// Pseudo-rapidity
          double eta() const;
          /// Rapidity
          double rapidity() const;
          void truncate( double tolerance = 1.e-10 );
          /// Rotate the transverse components by an angle phi (and reflect the y coordinate)
          Momentum& rotatePhi( double phi, double sign );
          /// Rotate the particle's momentum by a polar/azimuthal angle
          Momentum& rotateThetaPhi( double theta_, double phi_ );
          /// Apply a \f$ z\rightarrow -z\f$ transformation
          inline Momentum& mirrorZ() { pz_ = -pz_; return *this; }
        private:
          /// Compute the 3-momentum's norm
          void computeP();
          /// Momentum along the \f$x\f$-axis
          double px_;
          /// Momentum along the \f$y\f$-axis
          double py_;
          /// Momentum along the \f$z\f$-axis
          double pz_;
          /// 3-momentum's norm (in GeV/c)
          double p_;
          /// Energy (in GeV)
          double energy_;
      };
      /// Human-readable format for a particle's PDG code
      friend std::ostream& operator<<( std::ostream& os, const PDG& pc );
      /// Human-readable format for a particle's role in the event
      friend std::ostream& operator<<( std::ostream& os, const Particle::Role& rl );
      /// Compute the 4-vector sum of two 4-momenta
      friend Particle::Momentum operator+( const Particle::Momentum& mom1, const Particle::Momentum& mom2 );
      /// Compute the 4-vector difference of two 4-momenta
      friend Particle::Momentum operator-( const Particle::Momentum& mom1, const Particle::Momentum& mom2 );
      /// Compute the inverse per-coordinate 4-vector
      friend Particle::Momentum operator-( const Particle::Momentum& mom );
      /// Scalar product of two 3-momenta
      friend double operator*( const Particle::Momentum& mom1, const Particle::Momentum& mom2 );
      /// Multiply all components of a 4-momentum by a scalar
      friend Particle::Momentum operator*( const Particle::Momentum& mom, double c );
      /// Multiply all components of a 4-momentum by a scalar
      friend Particle::Momentum operator*( double c, const Particle::Momentum& mom );

      //----- static getters

      /// Convert a polar angle to a pseudo-rapidity
      static double thetaToEta( double theta );
      /// Convert a pseudo-rapidity to a polar angle
      static double etaToTheta( double eta );
      /// Convert a pseudo-rapidity to a rapidity
      static double etaToY( double eta_, double m_, double pt_ );

      Particle();
      /// Build using the role of the particle in the process and its PDG id
      /// \param[in] pdgId PDG identifier
      /// \param[in] role Role of the particle in the process
      /// \param[in] st Current status
      Particle( Role role, PDG pdgId, Status st = Status::Undefined );
      /// Copy constructor
      Particle( const Particle& );
      inline ~Particle() {}
      /// Comparison operator (from unique identifier)
      bool operator<( const Particle& rhs ) const;
      /// Comparison operator (from their reference's unique identifier)
      //bool operator<( Particle *rhs ) const { return ( id < rhs->id ); }

      // --- general particle properties

      /// Unique identifier (in a Event object context)
      int id() const { return id_; }
      //void setId( int id ) { id_ = id; }
      /// Set the particle unique identifier in an event
      void setId( int id ) { id_ = id; }
      /// Electric charge (given as a float number, for the quarks and bound states)
      float charge() const { return charge_sign_ * particleproperties::charge( pdg_id_ ); }
      /// Set the electric charge sign (+-1 for charged or 0 for neutral particles)
      void setChargeSign( int sign ) { charge_sign_ = sign; }
      /// Role in the considered process
      Role role() const { return role_; }
      /// Set the particle role in the process
      void setRole( const Role& role ) { role_ = role; }
      /**
       * Codes 1-10 correspond to currently existing partons/particles, and larger codes contain partons/particles which no longer exist, or other kinds of event information
       * \brief Particle status
       */
      Status status() const { return status_; }
      /// Set the particle decay/stability status
      void setStatus( Status status ) { status_ = status; }

      /// Set the PDG identifier (along with the particle's electric charge)
      /// \param[in] pdg PDG identifier
      /// \param[in] ch Electric charge (0, 1, or -1)
      void setPdgId( const PDG& pdg, short ch = 0 );
      /// Set the PDG identifier (along with the particle's electric charge)
      /// \param[in] pdg_id PDG identifier (incl. electric charge in e)
      void setPdgId( short pdg_id );
      /// Retrieve the objectified PDG identifier
      inline PDG pdgId() const { return pdg_id_; }
      /// Retrieve the integer value of the PDG identifier
      int integerPdgId() const;
      /// Particle's helicity
      float helicity() const { return helicity_; }
      /// Set the helicity of the particle
      void setHelicity( float heli ) { helicity_ = heli; }
      /// Particle mass in GeV/c\f$^2\f$
      /// \return Particle's mass
      inline double mass() const { return mass_; };
      /// Compute the particle mass
      /// \param[in] off_shell Allow the particle to be produced off-shell?
      /// \note This method ensures that the kinematics is properly set (the mass is set according to the energy and the momentum in priority)
      void computeMass( bool off_shell = false );
      /// Set the particle mass, in GeV/c\f$^2\f$
      /// \param m Mass in GeV/c\f$^2\f$
      /// \note This method ensures that the kinematics is properly set (the mass is set according to the energy and the momentum in priority)
      void setMass( double m = -1. );
      /// Particle squared mass, in GeV\f$^2\f$/c\f$^4\f$
      inline double mass2() const { return mass_*mass_; };
      /// Retrieve the momentum object associated with this particle
      inline Momentum& momentum() { return momentum_; }
      /// Retrieve the momentum object associated with this particle
      inline Momentum momentum() const { return momentum_; }
      /// Associate a momentum object to this particle
      void setMomentum( const Momentum& mom, bool offshell = false );
      /**
       * \brief Set the 3- or 4-momentum associated to the particle
       * \param[in] px Momentum along the \f$x\f$-axis, in GeV/c
       * \param[in] py Momentum along the \f$y\f$-axis, in GeV/c
       * \param[in] pz Momentum along the \f$z\f$-axis, in GeV/c
       * \param[in] e Energy, in GeV
       */
      void setMomentum( double px, double py, double pz, double e = -1. );
      /// Set the 4-momentum associated to the particle
      /// \param[in] p 4-momentum
      inline void setMomentum( double p[4] ) { setMomentum( p[0], p[1], p[2], p[3] ); }
      /// Set the particle's energy
      /// \param[in] e Energy, in GeV
      void setEnergy( double e = -1. );
      /// Get the particle's energy, in GeV
      double energy() const;
      /// Get the particle's squared energy, in GeV\f$^2\f$
      inline double energy2() const { return energy()*energy(); };
      /// Is this particle a valid particle which can be used for kinematic computations?
      bool valid();

      // --- particle relations

      /// Is this particle a primary particle?
      inline bool primary() const { return mothers_.empty(); }
      /// Set the mother particle
      /// \param[in] part A Particle object containing all the information on the mother particle
      void addMother( Particle& part );
      /// Get the unique identifier to the mother particle from which this particle arises
      /// \return An integer representing the unique identifier to the mother of this particle in the event
      inline ParticlesIds mothers() const { return mothers_; }
      /**
       * \brief Add a decay product
       * \param[in] part The Particle object in which this particle will desintegrate or convert
       * \return A boolean stating if the particle has been added to the daughters list or if it was already present before
       */
      void addDaughter( Particle& part );
      /// Gets the number of daughter particles
      inline unsigned int numDaughters() const { return daughters_.size(); };
      /// Get an identifiers list all daughter particles
      /// \return An integer vector containing all the daughters' unique identifier in the event
      inline ParticlesIds daughters() const { return daughters_; }

      // --- global particle information extraction

      /// Dump all the information on this particle into the standard output stream
      void dump() const;

    private:
      /// Unique identifier in an event
      int id_;
      /// Electric charge (+-1 or 0)
      short charge_sign_;
      /// Momentum properties handler
      Momentum momentum_;
      /// Mass, in GeV/c\f$^2\f$
      double mass_;
      /// Helicity
      float helicity_;
      /// Role in the process
      Role role_;
      /// Decay/stability status
      Status status_;
      /// List of mother particles
      ParticlesIds mothers_;
      /// List of daughter particles
      ParticlesIds daughters_;
      /// PDG id
      PDG pdg_id_;
  };

  /// Compute the centre of mass energy of two particles (incoming or outgoing states)
  double CMEnergy( const Particle& p1, const Particle& p2 );
  /// Compute the centre of mass energy of two particles (incoming or outgoing states)
  double CMEnergy( const Particle::Momentum& m1, const Particle::Momentum& m2 );

  //bool operator<( const Particle& a, const Particle& b ) { return a.id<b.id; }

  // --- particle containers

  /// List of Particle objects
  typedef std::vector<Particle> Particles;
  /// List of particles' roles
  typedef std::vector<Particle::Role> ParticleRoles;
  /// Map between a particle's role and its associated Particle object
  typedef std::unordered_map<Particle::Role,Particles> ParticlesMap;
}

#endif
