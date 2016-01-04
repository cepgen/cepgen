#ifndef Particle_h
#define Particle_h

#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <set>
#include <map>
#include <vector>
#include <algorithm>

#include "utils.h"
//#include "LHEutils.h"

typedef std::set<int> ParticlesIds;

/**
 * Kinematic information for one particle
 * @brief Kinematics of one particle
 */
class Particle {
  public:
    /** Unique identifier for a particle type. From @cite Beringer:1900zz :
     * _The Monte Carlo particle numbering scheme [...] is intended to facilitate interfacing between event generators, detector simulators, and analysis packages used in particle physics._
     * \brief PDG ids of all known particles
     */
    enum ParticleCode {
      invalidParticle = 0,
      dQuark = 1,
      uQuark = 2,
      Electron = 11,
      ElectronNeutrino = 12,
      Muon = 13,
      MuonNeutrino = 14,
      Tau = 15,
      TauNeutrino = 16,
      Gluon = 21,
      Photon = 22,
      PiPlus = 211,
      PiZero = 111,
      Rho770_0 = 113,
      Omega782 = 223,
      JPsi= 443,
      Phi1680 = 100333,
      Upsilon1S = 553,
      Upsilon2S = 100553,
      Upsilon3S = 200553,
      ud0Diquark = 2101,
      ud1Diquark = 2103,
      uu1Diquark = 2203,
      Proton = 2212,
      Neutron = 2112,
      Pomeron = 990,
      Reggeon = 110
    };
    enum Status {
      PrimordialIncoming = -9,
      Incoming = -1,
      Undecayed = -3,
      sPropagator = -2,
      Undefined = 0,
      FinalState = 1,
      Resonance = 2,
      DebugResonance = 3,
      PythiaHIncoming = 21,
      HerwigFragment = 193 //FIXME
    };
    /// Role of the particle in the process
    enum Role {
      UnknownRole = -1,
      IncomingBeam1 = 1,
      IncomingBeam2 = 2,
      Parton1 = 41,
      Parton2 = 42,
      Parton3 = 43,
      CentralSystem = 4,
      OutgoingBeam1 = 3,
      OutgoingBeam2 = 5,
      CentralParticle1 = 6,
      CentralParticle2 = 7
    };
    class Momentum {
      public:
        /// Build a 4-momentum at rest with an invalid energy (no mass information known)
        inline Momentum() : fPx(0.), fPy(0.), fPz(0.), fP(0.), fE(-1.) {;}
        /// Build a 4-momentum using its 3-momentum coordinates and its energy
        inline Momentum(double x_, double y_, double z_, double t_=-1.) :
          fPx(x_), fPy(y_), fPz(z_), fE(t_) { ComputeP(); }
        inline ~Momentum() {;}
        
        // --- static definitions
        
        /// Build a 3-momentum from its three pseudo-cylindric coordinates
        static inline Momentum FromPtEtaPhi(double pt, double eta, double phi, double e=-1.) {
          double px = pt*cos(phi), py = pt*sin(phi), pz = pt*sinh(eta);
          return Momentum(px, py, pz, e);
        }
        static inline Momentum FromPThetaPhi(double p, double theta, double phi, double e=-1.) {
          double px = p*sin(theta)*cos(phi), py = p*sin(theta)*sin(phi), pz = p*cos(theta);
          //std::cout << "sin(theta)=" << sin(theta) << ", sin(phi)=" << sin(phi) << std::endl;
          return Momentum(px, py, pz, e);
        }
        /// Build a 4-momentum from its four momentum and energy coordinates
        static inline Momentum FromPxPyPzE(double px, double py, double pz, double e) {
          return Momentum(px, py, pz, e);
        }
        
        // --- vector and scalar operators

        /// Add a 4-momentum through a 4-vector sum
        Momentum& operator+=(const Momentum&);
        /// Subtract a 4-momentum through a 4-vector sum
        Momentum& operator-=(const Momentum&);
        /// Scalar product of the 3-momentum with another 3-momentum
        double operator*=(const Momentum&);
        /// Multiply all 4-momentum coordinates by a scalar
        Momentum& operator*=(double c);
        /// Compute the 4-vector sum of two 4-momenta
        inline Momentum& operator+(const Momentum& mom_) { return *this += mom_; }
        /// Compute the 4-vector difference of two 4-momenta
        inline Momentum& operator-(const Momentum& mom_) { return *this -= mom_; }
        /// Scalar product of two 3-momenta
        inline double operator*(const Momentum& mom_) { return *this *= mom_; }
        /// Multiply all components of a 4-momentum by a scalar
        inline Momentum& operator*(double c) { return *this *= c; }

        void BetaGammaBoost(double gamma, double betagamma);

        // --- setters and getters
        
        /// Set all the components of the 4-momentum (in GeV)
        inline bool SetP(double px_,double py_,double pz_, double e_) {
          SetP(px_, py_, pz_); SetE(e_);
          return true;
        }
        /// Set all the components of the 3-momentum (in GeV)
        inline void SetP(double px_,double py_,double pz_) {
          fPx = px_;
          fPy = py_;
          fPz = pz_;
          ComputeP();
        }
        /// Set an individual component of the 4-momentum (in GeV)
        inline void SetP(unsigned int i, double p_) {
          switch (i) {
            case 0: fPx = p_; break;
            case 1: fPy = p_; break;
            case 2: fPz = p_; break;
            case 3: fE = p_; break;
            default: return;
          }
          ComputeP();
        }
        /// Set the energy (in GeV)
        inline void SetE(double e_) { fE = e_; }
        /// Get one component of the 4-momentum (in GeV)
        inline double P(unsigned int i) const {
          switch (i) {
            case 0: return fPx;
            case 1: return fPy;
            case 2: return fPz;
            case 3: return fE;
            default: return -1.;
          }
        }
        /// Get the momentum along the \f$x\f$-axis (in GeV)
        inline double Px() const { return fPx; }
        /// Get the momentum along the \f$y\f$-axis (in GeV)
        inline double Py() const { return fPy; }
        /// Get the longitudinal momentum (in GeV)
        inline double Pz() const { return fPz; }
        /// Get the transverse momentum (in GeV)
        inline double Pt() const { return sqrt(pow(Px(),2)+pow(Py(),2)); }
        /// Get the 3-momentum norm (in GeV)
        inline double P() const { return fP; }
        /// Get the squared 3-momentum norm (in \f$\text{GeV}^\text{2}\f$)
        inline double P2() const { return pow(fP,2); }
        /// Get the energy (in GeV)
        inline double E() const { return fE; }
        /// Get the particle's mass (in GeV) as computed from its energy and momentum
        inline double M() const { return sqrt(pow(E(),2)-P2()); }
        inline double Theta() const { return atan2(Pt(), Pz()); }
        /// Get the azimutal angle (angle in the transverse plane)
        inline double Phi() const { return atan2(Py(), Px()); }
        /// Get the pseudo-rapidity
        inline double Eta() const {
          return (Pt()!=0.)
            ? log((P()+fabs(Pz()))/Pt())*(Pz()/fabs(Pz()))
            : 9999.*(Pz()/fabs(Pz()));
        };
        /// Get the rapidity
        inline double Rapidity() const { return (E()<0.) ? 999. : log((E()+Pz())/(E()-Pz()))/2.; };
        /// Rotate the transverse components by an angle phi (and reflect the y coordinate)
        inline void RotatePhi(double phi, double rany) {
          double px = fPx*cos(phi)+fPy*sin(phi)*rany,
                 py =-fPx*sin(phi)+fPy*cos(phi)*rany;
          fPx = px;
          fPy = py;
        }
      private:
        /// Compute the 3-momentum's norm
        inline void ComputeP() {
          fP = 0.;
          for (unsigned int i=0; i<3; i++) { fP += pow(P(i),2); }
          fP = sqrt(fP);
        }
        /// Momentum along the \f$x\f$-axis
        double fPx;
        /// Momentum along the \f$y\f$-axis
        double fPy;
        /// Momentum along the \f$z\f$-axis
        double fPz;
        /// 3-momentum's norm (in GeV/c)
        double fP;
        /// Energy (in GeV)
        double fE;
    };
    friend std::ostream& operator<<(std::ostream& os, const Particle::ParticleCode& pc);
    /**
     * Gets the mass in GeV/c**2 of a particle given its PDG identifier
     * @brief Gets the mass of a particle
     * @param pdgId_ ParticleCode (PDG ID)
     * @return Mass of the particle in \f$\text{GeV}/c^2\f$
     */
    static double GetMassFromPDGId(Particle::ParticleCode pdgId_);
    /**
     * Gets the total decay width for one particle to be decayed
     * @param[in] ParticleCode (PDG ID)
     * @return Decay width in GeV
     */
    static double GetWidthFromPDGId(Particle::ParticleCode pdgId_);
    
    Particle();
    /// Build using the role of the particle in the process and its PDG id
    Particle(Role role_,ParticleCode pdgId_=Particle::invalidParticle);
    inline ~Particle() {;}
    /// Assignment operator
    Particle& operator=(const Particle&);
    /// Comparison operator (from unique identifier)
    inline bool operator<(const Particle& rhs) { return id<rhs.id; };
    /// Comparison operator (from their reference's unique identifier)
    inline bool operator<(const Particle *rhs) { std::cout << id << "\t" << rhs->id << "\t" << (id<rhs->id) << std::endl; return id<rhs->id; };
    void LorentzBoost(double m_, const Momentum& mom_);
    /// Lorentz boost (from ROOT)
    double* LorentzBoost(const Momentum& mom_);
    
    // --- general particle properties
    
    /// Unique identifier (in a Event object context)
    int id;
    /// Electric charge (given as a float number, for the quarks and bound states)
    float charge;
    /// Human-readable name
    std::string name;
    /// Role in the considered process
    Role role;
    /**
     * Codes 1-10 correspond to currently existing partons/particles, and larger codes contain partons/particles which no longer exist, or other kinds of event information
     * @brief Particle status
     */
    Status status;

    inline void SetPDGId(ParticleCode pdg, float ch=-999.) {
      fPDGid = pdg;
      if (ch==-999.) charge = pdg/abs(pdg);
      else charge = ch;
    }
    inline ParticleCode GetPDGId() const { return fPDGid; }
    inline int GetIntPDGId() const {
      int pdg = static_cast<int>(fPDGid);
      if (pdg>10 and pdg<16 and pdg%2!=0) return static_cast<int>(-charge)*pdg;
      else return pdg;
    }
    /// Particle's helicity
    /// @fixme Float??
    float helicity;
    /**
     * Gets the particle's mass in \f$\text{GeV}/c^{2}\f$.
     * @brief Gets the particle's mass
     * @return The particle's mass
     */
    inline double M() const { return fMass; };
    /**
     * Set the mass of the particle in \f$\text{GeV}/c^{2}\f$ according to a value given as an argument. This method ensures that the kinematics is properly set (the mass is set according to the energy and the momentum in priority)
     * @param m_ The mass in \f$\text{GeV}/c^{2}\f$ to set
     * @brief Set the particle's mass in \f$\text{GeV}/c^{2}\f$
     * @return A boolean stating whether or not the mass was correctly set
     */
    bool SetM(double m_=-1.);
    /// Get the particle's squared mass (in \f$\text{GeV}^\text{2}\f$)
    inline double M2() const { return std::pow(fMass,2); };
    inline Momentum GetMomentum() const { return fMomentum; }
    inline bool SetMomentum(const Momentum& mom) {
      fMomentum = mom;
      if (fMass<0.) { SetM(); }
      double e = sqrt(fMomentum.P2()+pow(fMass,2));
      if (mom.E()<0.) {
        fMomentum.SetE(e);
        return true;
      }
      if (fabs(e-fMomentum.E())<1.e-6) { // less than 1 eV difference
        return true;
      }
      if (fabs(e-mom.E())<1.e-6) { // less than 1 eV difference
        return true;
      }
      if (role!=Parton1 and role!=Parton2) {
        Error(Form("Energy difference for particle %d (computed-set): %.5f", (int)role, e-fMomentum.E()));
      }
      fMomentum.SetE(e);
      return false;
    }
    /**
     * @brief Set the 3-momentum associated to the particle
     * @param[in] px_ Momentum along the \f$x\f$-axis, in \f$\text{GeV}/c\f$
     * @param[in] py_ Momentum along the \f$y\f$-axis, in \f$\text{GeV}/c\f$
     * @param[in] pz_ Momentum along the \f$z\f$-axis, in \f$\text{GeV}/c\f$
     * @return A boolean stating the validity of this particle (according to its 4-momentum norm)
     */
    inline bool SetMomentum(double px_,double py_,double pz_) {
      fMomentum.SetP(px_, py_, pz_); SetE();
      return true;
    };
    /**
     * @brief Set the 4-momentum associated to the particle
     * @param[in] px_ Momentum along the \f$x\f$-axis, in \f$\text{GeV}/c\f$
     * @param[in] py_ Momentum along the \f$y\f$-axis, in \f$\text{GeV}/c\f$
     * @param[in] pz_ Momentum along the \f$z\f$-axis, in \f$\text{GeV}/c\f$
     * @param[in] e_ Energy, in GeV
     * @return A boolean stating the validity of the particle's kinematics
     */
    inline bool SetMomentum(double px_,double py_,double pz_,double e_) {
      SetMomentum(px_, py_, pz_);
      if (fabs(e_-fMomentum.E())>1.e-6) { // more than 1 eV difference
        Error(Form("Energy difference: %.5f", e_-fMomentum.E()));
        return false;
      }
      return true;
    };
    /**
     * @brief Set the 4-momentum associated to the particle
     * @param[in] p_ 4-momentum
     * @return A boolean stating the validity of the particle's kinematics
     */
    inline bool SetMomentum(double p_[4]) { return SetMomentum(p_[0], p_[1], p_[2], p_[3]); };
    /**
     * @brief Set the particle's energy
     * @param[in] e_ Energy, in GeV
     */
    inline void SetE(double e_=-1.) {
      if (e_<0. and fMass>=0.) e_ = sqrt(M2()+fMomentum.P2());
      fMomentum.SetE(e_);
    }
    /// Get the particle's energy (in GeV)
    inline double E() const {
      return (fMomentum.E()<0.) ? std::sqrt(M2()+fMomentum.P2()) : fMomentum.E();
    };
    /// Get the particle's squared energy (in \f$\text{GeV}^\text{2}\f$)
    inline double E2() const { return std::pow(E(), 2); };
    void RotateThetaPhi(double theta_, double phi_);
    /// Is this particle a valid particle which can be used for kinematic computations ?
    bool Valid();
    
    // --- particle relations
    
    /// Is this particle a primary particle ?
    inline bool Primary() const { return fIsPrimary; }
    /**
     * @brief Set the mother particle
     * @param[in] part_ A Particle object containing all the information on the mother particle
     */
    void SetMother(Particle* part_);
    /**
     * @brief Gets the unique identifier to the mother particle from which this particle arises
     * @return An integer representing the unique identifier to the mother of this particle in the event
     */
    inline ParticlesIds GetMothersIds() const { return fMothers; }
    /**
     * @brief Add a decay product
     * @param[in] part_ The Particle object in which this particle will desintegrate or convert
     * @return A boolean stating if the particle has been added to the daughters list or if it was already present before
     */
    bool AddDaughter(Particle* part_);
    /// Gets the number of daughter particles
    inline unsigned int NumDaughters() { return fDaughters.size(); };
    /**
     * @brief Get an identifiers list all daughter particles
     * @return An integer vector containing all the daughters' unique identifier in the event
     */
    std::vector<int> GetDaughters();
    
    // --- global particle information extraction
    
    /**
     * Returns a string containing all the particle's kinematics as expressed in the Les Houches format
     * @param[in] revert_ Is the event symmetric ? If set to true, the third component of the momentum is reverted.
     * @return The LHE line associated to the particle, and containing the particle's history (mother/daughters), its kinematics, and its status
     */
    std::string GetLHEline(bool revert_=false);
    /// Dump all the information on this particle into the standard output stream
    void Dump();    
    void PDF2PDG();
    
    // --- other methods
    
    /**
     * Hadronise the particle with a generic hadroniser, and builds the shower (list of Particle objects) embedded in this object
     * @param[in] algo_ Algorithm in use to hadronise the particle
     * @brief Hadronises the particle using Pythia
     * @return A boolean stating whether or not the particle has been hadronised
     */
    bool Hadronise(std::string algo_);
  private:
    Momentum fMomentum;
    /// Mass in \f$\text{GeV}/c^2\f$
    double fMass;
    /// List of mother particles
    ParticlesIds fMothers;
    /// List of daughter particles
    ParticlesIds fDaughters;
    /// PDG id
    ParticleCode fPDGid;
    /// Is the particle a primary particle ?
    bool fIsPrimary;
    double __tmp3[3], __tmp4[4];
};

inline bool compareParticle(Particle a, Particle b) { return a.id<b.id; }
inline bool compareParticlePtrs(Particle* a, Particle* b) { return a->id<b->id; }

/// Compute the centre of mass energy of two particles (incoming or outgoing states)
inline static double CMEnergy(const Particle& p1, const Particle& p2) {
  if (p1.M()*p2.M()<0.) return 0.;
  if (p1.E()*p2.E()<0.) return 0.;
  return sqrt(p1.M2()+p2.M2()+2.*p1.E()*p2.E()-2.*(p1.GetMomentum()*p2.GetMomentum()));
}

// --- particle containers

typedef std::vector<Particle> Particles;
typedef std::vector<Particle*> ParticlesRef;
typedef std::vector<Particle::Role> ParticleRoles;
typedef std::multimap<Particle::Role,Particle> ParticlesMap;

#endif
