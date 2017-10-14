#ifndef CepGen_Event_Particle_h
#define CepGen_Event_Particle_h

#include <set>
#include <map>
#include <vector>

#include "CepGen/Core/utils.h"

namespace CepGen
{

  /// A set of integer-type particle identifiers
  typedef std::set<int> ParticlesIds;

  /// Kinematic information for one particle
  class Particle {
    public:
      /** Unique identifier for a particle type. From \cite Beringer:1900zz :
       * `The Monte Carlo particle numbering scheme [...] is intended to facilitate interfacing between event generators, detector simulators, and analysis packages used in particle physics.`
       * \brief PDG ids of all known particles
       */
      enum ParticleCode {
        invalidParticle = 0,
        //--- fundamental particles
        dQuark = 1, uQuark = 2,
        Electron = 11, ElectronNeutrino = 12,
        Muon = 13, MuonNeutrino = 14,
        Tau = 15, TauNeutrino = 16,
        Gluon = 21, Photon = 22, Z = 23, W = 24,
        //--- composite particles
        PiPlus = 211, PiZero = 111,
        Rho770_0 = 113, Rho1450_0 = 100113, Rho1700_0 = 30113,
        Eta = 221, Omega782 = 223,
        h1380_1 = 10333,
        JPsi= 443,
        Phi1680 = 100333,
        Upsilon1S = 553, Upsilon2S = 100553, Upsilon3S = 200553,
        ud0Diquark = 2101, ud1Diquark = 2103, uu1Diquark = 2203,
        Proton = 2212, Neutron = 2112,
        Pomeron = 990, Reggeon = 110
      };
      /// Internal status code for a particle
      enum Status {
        PrimordialIncoming = -9,
        Undecayed = -3,
        Propagator = -2,
        Incoming = -1,
        Undefined = 0,
        FinalState = 1,
        Resonance = 2,
        DebugResonance = 3,
        PythiaHIncoming = 21
      };
      /// Role of the particle in the process
      enum Role {
        UnknownRole = -1,
        IncomingBeam1 = 1, IncomingBeam2 = 2,
        OutgoingBeam1 = 3, OutgoingBeam2 = 5,
        CentralSystem = 6,
        Intermediate = 4,
        Parton1 = 41, Parton2 = 42, Parton3 = 43
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
          Momentum( double x_, double y_, double z_, double t_=-1. );
          inline ~Momentum() {}

          // --- static definitions

          /// Build a 3-momentum from its three pseudo-cylindric coordinates
          static Momentum fromPtEtaPhi( double pt, double eta, double phi, double e=-1. );
          /// Build a 4-momentum from its scalar momentum, and its polar and azimuthal angles
          static Momentum fromPThetaPhi( double p, double theta, double phi, double e=-1. );
          /// Build a 4-momentum from its four momentum and energy coordinates
          static Momentum fromPxPyPzE( double px, double py, double pz, double e );

          // --- vector and scalar operators

          /// Scalar product of the 3-momentum with another 3-momentum
          double threeProduct( const Momentum& ) const;
          /// Scalar product of the 4-momentum with another 4-momentum
          double fourProduct( const Momentum& ) const;
          /// Add a 4-momentum through a 4-vector sum
          Momentum& operator+=( const Momentum& );
          /// Subtract a 4-momentum through a 4-vector sum
          Momentum& operator-=( const Momentum& );
          /// Scalar product of the 3-momentum with another 3-momentum
          double operator*=( const Momentum& );
          /// Multiply all 4-momentum coordinates by a scalar
          Momentum& operator*=( double c );
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
          /// Set an individual component of the 4-momentum (in GeV)
          void setP( unsigned int i, double p );
          /// Set the energy (in GeV)
          inline void setEnergy( double e ) { energy_ = e; }
          /// Compute the energy from the mass
          inline void setMass( double m ) { energy_ = sqrt( p2()+m*m ); }
          /// Compute the energy from the mass
          inline void setMass2( double m2 ) { energy_ = sqrt( p2()+m2 ); }
          /// Get one component of the 4-momentum (in GeV)
          double operator[]( const unsigned int i ) const;
          /// Momentum along the \f$x\f$-axis (in GeV)
          inline double px() const { return px_; }
          /// Momentum along the \f$y\f$-axis (in GeV)
          inline double py() const { return py_; }
          /// Longitudinal momentum (in GeV)
          inline double pz() const { return pz_; }
          /// Transverse momentum (in GeV)
          inline double pt() const { return sqrt( pt2() ); }
          /// Squared transverse momentum (in GeV\f$^\textrm{2}\f$)
          inline double pt2() const { return ( px()*px()+py()*py() ); }
          /// 4-vector of double precision floats (in GeV)
          const std::vector<double> pVector() const;
          /// 3-momentum norm (in GeV)
          inline double p() const { return p_; }
          /// Squared 3-momentum norm (in GeV\f$^\textrm{2}\f$)
          inline double p2() const { return p_*p_; }
          /// Energy (in GeV)
          inline double energy() const { return energy_; }
          /// Squared energy (in GeV^2)
          inline double energy2() const { return energy_*energy_; }
          /// Squared mass (in GeV^2) as computed from its energy and momentum
          inline double mass2() const { return energy2()-p2(); }
          /// Mass (in GeV) as computed from its energy and momentum
          /// \note Returns \f$-\sqrt{|E^2-\mathbf{p}^2|}<0\f$ if \f$\mathbf{p}^2>E^2\f$
          double mass() const;
          /// Polar angle (angle with respect to the longitudinal direction)
          inline double theta() const { return atan2( pt(), pz() ); }
          /// Azimutal angle (angle in the transverse plane)
          inline double phi() const { return atan2( py(), px() ); }
          /// Pseudo-rapidity
          double eta() const;
          /// Rapidity
          double rapidity() const;
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
      friend std::ostream& operator<<( std::ostream& os, const Particle::ParticleCode& pc );
      /// Human-readable format for a particle's role in the event
      friend std::ostream& operator<<( std::ostream& os, const Particle::Role& rl );
      /// Compute the 4-vector sum of two 4-momenta
      friend Particle::Momentum operator+( const Particle::Momentum& mom1, const Particle::Momentum& mom2 );
      /// Compute the 4-vector difference of two 4-momenta
      friend Particle::Momentum operator-( const Particle::Momentum& mom1, const Particle::Momentum& mom2 );
      /// Scalar product of two 3-momenta
      friend double operator*( const Particle::Momentum& mom1, const Particle::Momentum& mom2 );
      /// Multiply all components of a 4-momentum by a scalar
      friend Particle::Momentum operator*( const Particle::Momentum& mom, double c );
      /// Multiply all components of a 4-momentum by a scalar
      friend Particle::Momentum operator*( double c, const Particle::Momentum& mom );

      //----- static getters

      /**
       * \brief Gets the mass of a particle
       * \param pdgId ParticleCode (PDG ID)
       * \return Mass of the particle in \f$\textrm{GeV}/c^2\f$
       */
      static double massFromPDGId( const Particle::ParticleCode& pdgId );
      /**
       * \brief Gets the electric charge of a particle
       * \param pdgId ParticleCode (PDG ID)
       * \return Charge of the particle in \f$e\f$
       */
      static double chargeFromPDGId( const Particle::ParticleCode& pdgId );
      /**
       * \brief Total decay width of one unstable particle
       * \param[in] pdgId ParticleCode (PDG ID)
       * \return Decay width in GeV
       */
      static double widthFromPDGId( const Particle::ParticleCode& pdgId );

      /// Convert a polar angle to a pseudo-rapidity
      static inline double thetaToEta( double theta_ ) { return -log( tan( theta_/180.*M_PI/2. ) ); }
      /// Convert a pseudo-rapidity to a polar angle
      static inline double etaToTheta( double eta_ ) { return 2.*atan( exp( -eta_ ) )*180. / M_PI; }
      /// Convert a pseudo-rapidity to a rapidity
      static double etaToY( double eta_, double m_, double pt_ );



      Particle();
      /// Build using the role of the particle in the process and its PDG id
      /// \param[in] pdgId ParticleCode (PDG ID)
      /// \param[in] role Role of the particle in the process
      /// \param[in] st Current status
      Particle( Role role, ParticleCode pdgId=invalidParticle, Status st=Undefined );
      /// Copy constructor
      Particle( const Particle& );
      inline ~Particle() {}
      /// Comparison operator (from unique identifier)
      bool operator<( const Particle& rhs ) const;
      /// Comparison operator (from their reference's unique identifier)
      //bool operator<( Particle *rhs ) const { return ( id < rhs->id ); }
      Particle& lorentzBoost( double m_, const Momentum& mom_ );
      /// Lorentz boost (shamelessly stolen from ROOT)
      std::vector<double> lorentzBoost( const Momentum& mom_ );

      // --- general particle properties

      /// Unique identifier (in a Event object context)
      int id() const { return id_; }
      //void setId( int id ) { id_ = id; }
      /// Set the particle unique identifier in an event
      void setId( int id ) { id_ = id; }
      /// Electric charge (given as a float number, for the quarks and bound states)
      float charge() const { return charge_sign_ * chargeFromPDGId( pdg_id_ ); }
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
      void setStatus( const Status& status ) { status_ = status; }

      /// Set the PDG identifier (along with the particle's electric charge)
      /// \param[in] pdg ParticleCode (PDG ID)
      /// \param[in] ch Electric charge (0, 1, or -1)
      void setPdgId( const ParticleCode& pdg, short ch=0 );
      /// Retrieve the objectified PDG identifier
      inline ParticleCode pdgId() const { return pdg_id_; }
      /// Retrieve the integer value of the PDG identifier
      int integerPdgId() const;
      /// Particle's helicity
      float helicity() const { return helicity_; }
      /// Set the helicity of the particle
      void setHelicity( float heli ) { helicity_ = heli; }
      /**
       * Gets the particle's mass in \f$\textrm{GeV}/c^{2}\f$.
       * \brief Gets the particle's mass
       * \return The particle's mass
       */
      inline double mass() const { return mass_; };
      /**
       * Set the mass of the particle in \f$\textrm{GeV}/c^{2}\f$ while ensuring that the kinematics is properly set (the mass is set according to the energy and the momentum in priority)
       * \brief Compute the particle's mass in \f$\textrm{GeV}/c^{2}\f$
       */
      void computeMass( bool off_shell=false );
      /**
       * Set the mass of the particle in \f$\textrm{GeV}/c^{2}\f$ according to a value given as an argument. This method ensures that the kinematics is properly set (the mass is set according to the energy and the momentum in priority)
       * \param m The mass in \f$\textrm{GeV}/c^{2}\f$ to set
       * \brief Set the particle's mass in \f$\textrm{GeV}/c^{2}\f$
       */
      void setMass( double m=-1. );
      /// Get the particle's squared mass (in \f$\textrm{GeV}^\textrm{2}\f$)
      inline double mass2() const { return mass_*mass_; };
      /// Retrieve the momentum object associated with this particle
      inline Momentum& momentum() { return momentum_; }
      /// Retrieve the momentum object associated with this particle
      inline Momentum momentum() const { return momentum_; }
      /// Associate a momentum object to this particle
      void setMomentum( const Momentum& mom, bool offshell=false );
      /**
       * \brief Set the 3-momentum associated to the particle
       * \param[in] px Momentum along the \f$x\f$-axis, in \f$\textrm{GeV}/c\f$
       * \param[in] py Momentum along the \f$y\f$-axis, in \f$\textrm{GeV}/c\f$
       * \param[in] pz Momentum along the \f$z\f$-axis, in \f$\textrm{GeV}/c\f$
       */
      void setMomentum( double px, double py, double pz );
      /**
       * \brief Set the 4-momentum associated to the particle
       * \param[in] px Momentum along the \f$x\f$-axis, in \f$\textrm{GeV}/c\f$
       * \param[in] py Momentum along the \f$y\f$-axis, in \f$\textrm{GeV}/c\f$
       * \param[in] pz Momentum along the \f$z\f$-axis, in \f$\textrm{GeV}/c\f$
       * \param[in] e Energy, in GeV
       */
      void setMomentum( double px, double py, double pz, double e );
      /**
       * \brief Set the 4-momentum associated to the particle
       * \param[in] p 4-momentum
       */
      inline void setMomentum( double p[4] ) { setMomentum( p[0], p[1], p[2], p[3] ); }
      /**
       * \brief Set the particle's energy
       * \param[in] e Energy, in GeV
       */
      void setEnergy( double e=-1. );
      /// Get the particle's energy (in GeV)
      inline double energy() const {
        return ( momentum_.energy()<0. ) ? std::sqrt( mass2()+momentum_.p2() ) : momentum_.energy();
      };
      /// Get the particle's squared energy (in \f$\textrm{GeV}^\textrm{2}\f$)
      inline double energy2() const { return energy()*energy(); };
      /// Is this particle a valid particle which can be used for kinematic computations ?
      bool valid();

      // --- particle relations

      /// Is this particle a primary particle ?
      inline bool primary() const { return mothers_.empty(); }
      /**
       * \brief Set the mother particle
       * \param[in] part A Particle object containing all the information on the mother particle
       */
      void addMother( Particle& part );
      /**
       * \brief Gets the unique identifier to the mother particle from which this particle arises
       * \return An integer representing the unique identifier to the mother of this particle in the event
       */
      inline ParticlesIds mothers() const { return mothers_; }
      /**
       * \brief Add a decay product
       * \param[in] part The Particle object in which this particle will desintegrate or convert
       * \return A boolean stating if the particle has been added to the daughters list or if it was already present before
       */
      void addDaughter( Particle& part );
      /// Gets the number of daughter particles
      inline unsigned int numDaughters() const { return daughters_.size(); };
      /**
       * \brief Get an identifiers list all daughter particles
       * \return An integer vector containing all the daughters' unique identifier in the event
       */
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
      /// Mass in \f$\textrm{GeV}/c^2\f$
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
      ParticleCode pdg_id_;
  };

  /// Compute the centre of mass energy of two particles (incoming or outgoing states)
  inline static double CMEnergy( const Particle& p1, const Particle& p2 ) {
    if ( p1.mass()*p2.mass() < 0. ) return 0.;
    if ( p1.energy()*p2.energy() < 0. ) return 0.;
    return sqrt( p1.mass2()+p2.mass2() + 2.*p1.energy()*p2.energy() - 2.*( p1.momentum()*p2.momentum() ) );
  }

  /// Compute the centre of mass energy of two particles (incoming or outgoing states)
  inline static double CMEnergy( const Particle::Momentum& m1, const Particle::Momentum& m2 ) {
    if ( m1.mass()*m2.mass() < 0. ) return 0.;
    if ( m1.energy()*m2.energy() < 0. ) return 0.;
    return sqrt( m1.mass2()+m2.mass2() + 2.*m1.energy()*m2.energy() - 2.*( m1*m2 ) );
  }

  //bool operator<( const Particle& a, const Particle& b ) { return a.id<b.id; }

  // --- particle containers

  /// List of Particle objects
  typedef std::vector<Particle> Particles;
  /// List of particles' roles
  typedef std::vector<Particle::Role> ParticleRoles;
  /// Map between a particle's role and its associated Particle object
  typedef std::map<Particle::Role,Particles> ParticlesMap;
}

#endif
