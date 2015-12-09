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
    enum ParticleCode {
      invalidParticle = 0,
      dQuark = 1,
      uQuark = 2,
      Electron = 11,
      Muon = 13,
      Tau = 15,
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
      PythiaHIncoming = 21
    };
    class Momentum {
      public:
        inline Momentum() : fPx(0.), fPy(0.), fPz(0.), fE(-1.) {;}
        inline Momentum(double x_, double y_, double z_, double t_=-1.) :
          fPx(x_), fPy(y_), fPz(z_), fE(t_) {;}
        inline ~Momentum() {;}

        Momentum& operator+=(const Momentum&);
        Momentum& operator-=(const Momentum&);
        Momentum& operator-(const Momentum&);
        double operator*(const Momentum&); // scalar product

        inline bool SetP(double px_,double py_,double pz_, double e_) {
          SetP(px_, py_, pz_); SetE(e_);
          return true;
        }
        inline void SetP(double px_,double py_,double pz_) {
          fPx = px_; fPy = py_, fPz = pz_;
        }
        inline void SetP(unsigned int i, double p_) {
          switch (i) {
            case 0: fPx = p_; break;
            case 1: fPy = p_; break;
            case 2: fPz = p_; break;
            case 3: fE = p_; break;
            default: return;
          }
        }
        inline void SetE(double e_) { fE = e_; }
        inline double P(unsigned int i) const {
          switch (i) {
            case 0: return fPx;
            case 1: return fPy;
            case 2: return fPz;
            case 3: return fE;
            default: return -1.;
          }
        }
        inline double Px() const { return fPx; }
        inline double Py() const { return fPy; }
        inline double Pz() const { return fPz; }
        inline double Pt() const { return sqrt(pow(Px(),2)+Py()); }
        inline double P2() const { return pow(P(0),2)+pow(P(1),2)+pow(P(2),2); }
        inline double E() const { return fE; }
        inline double M() const { return sqrt(pow(E(),2)-P2()); }
        inline void RotatePhi(double phi, double rany) {
          double px = fPx*cos(phi)+fPy*sin(phi)*rany,
                 py =-fPx*sin(phi)+fPy*cos(phi)*rany;
          fPx = px;
          fPy = py;
        }
      private:
        double fPx;
        double fPy;
        double fPz;
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
    /**
     * @brief Object constructor (providing the role of the particle in the process, and its Particle Data Group identifier)
     */
    Particle(int role_,ParticleCode pdgId_=Particle::invalidParticle);
    ~Particle();
    /**
     * @brief Copies all the relevant quantities from one Particle object to another
     */
    Particle& operator=(const Particle&);
    /**
     * @brief Adds two particles' momenta to create a combined particle
     */
    Particle& operator+=(const Particle&);
    /**
     * @brief Substracts two particles' momenta to extract a particle's kinematics
     */
    Particle& operator-(const Particle&);
    /**
     * @brief Comparison operator to enable the sorting of particles in an event according to their unique identifier
     */
    inline bool operator<(const Particle& rhs) { std::cout << id << "\t" << rhs.id << "\t" << (id<rhs.id) << std::endl; return id<rhs.id; };
    /**
     * @brief Comparison operator to enable the sorting of Particle objects' pointers in an event according to their reference's unique identifier
     */
    inline bool operator<(const Particle *rhs) { std::cout << id << "\t" << rhs->id << "\t" << (id<rhs->id) << std::endl; return id<rhs->id; };
    /**
     * Returns a string containing all the particle's kinematics as expressed in the Les Houches format
     * @param[in] revert_ Is the event symmetric ? If set to true, the third component of the momentum is reverted.
     * @return The LHE line associated to the particle, and containing the particle's history (mother/daughters), its kinematics, and its status
     */
    std::string GetLHEline(bool revert_=false);
    /**
     * Dumps into the standard output stream all the available information on this particle
     * @brief Dumps all the information on this particle
     */
    void Dump();
    //double* LorentzBoost(double m_, double p_[4]);
    void LorentzBoost(double m_, const Momentum& mom_);
    /**
     * Lorentz boost from ROOT
     */
    double* LorentzBoost(const Momentum& mom_);
    /**
     * @brief Unique identifier of the particle (in a Event object context)
     */
    int id;
    /**
     * @brief The particle's electric charge (given as a float number, for the quarks and bound states)
     */
    float charge;
    /**
     * @brief Particle's name in a human-readable format
     */
    std::string name;
    /**
     * @brief Role in the considered process
     */
    int role;
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
    /**
     * Particle's helicity
     * @fixme Float??
     */
    float helicity;
    /**
     * @brief Momentum along the \f$x\f$-axis in \f$\text{GeV}/c\f$
     */
    inline double Px() const { return fMomentum.Px(); };
    /**
     * @brief Momentum along the \f$y\f$-axis in \f$\text{GeV}/c\f$
     */
    inline double Py() const { return fMomentum.Py(); };
    /**
     * @brief Momentum along the \f$z\f$-axis in \f$\text{GeV}/c\f$
     */
    inline double Pz() const { return fMomentum.Pz(); };
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
    /**
     * @brief Gets the particle's squared mass
     */
    inline double M2() const { return std::pow(fMass,2); };
    /**
     * Computes and returns \f$\eta\f$, the pseudo-rapidity of the particle
     * @brief Pseudo-rapidity
     * @return The pseudo-rapidity of the particle
     */
    inline double Eta() const {
      return (Pt()!=0.)
	      ? log((P()+fabs(Pz()))/Pt())*(Pz()/fabs(Pz()))
	      : 9999.*(Pz()/fabs(Pz()));
    };
    /**
     * Computes and returns \f$y\f$, the rapidity of the particle
     * @brief Rapidity
     * @return The rapidity of the particle
     */
    inline double Rapidity() { return (E()<0.) ? 999. : log((E()+Pz())/(E()-Pz()))/2.; };
    /**
     * Computes and returns \f$\phi\f$, the azimuthal angle of the particle in the transverse plane
     * @brief Azimuthal angle
     * @return The azimuthal angle of the particle
     */
    /*inline double Phi() const {
      return (Px()==0. && Py()==0. && Pz()==0.)
	      ? 0.
	      : atan2(P(), Pz());
    };*/
    inline void SetMomentum(const Momentum& mom) {
      fMomentum = mom;
      if (fMass>=0.) {
        double e = sqrt(fMomentum.P2()+pow(fMass,2));
        fMomentum.SetE(e);
      }
    }
    inline Momentum GetMomentum() const { return fMomentum; }
    /**
     * @brief Sets the 3-momentum associated to the particle
     * @param[in] px_ Momentum along the \f$x\f$-axis, in \f$\text{GeV}/c\f$
     * @param[in] py_ Momentum along the \f$y\f$-axis, in \f$\text{GeV}/c\f$
     * @param[in] pz_ Momentum along the \f$z\f$-axis, in \f$\text{GeV}/c\f$
     * @return A boolean stating the validity of this particle (according to its 4-momentum norm)
     */
    inline bool P(double px_,double py_,double pz_) {
      fMomentum.SetP(px_, py_, pz_);
      SetM(); SetE();
      return true;
    };
    /**
     * Sets the 4-momentum associated to the particle, and computes its (invariant) mass.
     * @brief Sets the 4-momentum associated to the particle
     * @param[in] px_ Momentum along the \f$x\f$-axis, in \f$\text{GeV}/c\f$
     * @param[in] py_ Momentum along the \f$y\f$-axis, in \f$\text{GeV}/c\f$
     * @param[in] pz_ Momentum along the \f$z\f$-axis, in \f$\text{GeV}/c\f$
     * @param[in] E_ Energy, in GeV
     * @return A boolean stating the validity of the particle's kinematics
     */
    inline bool P(double px_,double py_,double pz_,double E_) { double pp4[] = { px_, py_, pz_, E_ }; return P(pp4); };
    /**
     * @brief Sets the 4-momentum associated to the particle
     * @param[in] p_ 4-momentum
     * @return A boolean stating the validity of the particle's kinematics
     */
    inline bool P(double p_[4]) {
      fMomentum.SetP(p_[0], p_[1], p_[2], p_[3]);
      if (p_[3]<0.) return fMomentum.M()==fMass;
      else return true;
    };
    /**
     * Get one component of the particle's 4-momentum
     * @param[in] c_ The component to retrieve:
     * - 0-2: \f$\mathbf p = (p_x, p_y, p_z)\f$ (in \f$\text{GeV}/c\f$)
     * - 3: \f$E\f$ (in \f$\text{GeV}\f$)
     * @return The requested component of the energy-momentum for the particle
     */
    inline double P(int c_) const {
      if (c_>=0 and c_<4) return fMomentum.P(c_);
      else if (c_==4) return M();
      else return -999.;
    };
    /**
     * @brief Transverse momentum, in \f$\text{GeV}/c\f$
     */
    inline double Pt() const { return std::sqrt(std::pow(Px(),2)+std::pow(Py(),2)); };
    /**
     * @brief Norm of the 3-momentum, in \f$\text{GeV}/c\f$
     * @return The particle's 3-momentum norm as a double precision float
     */
    inline double P() const { return std::sqrt(std::pow(Pt(),2)+std::pow(Pz(),2)); };
    /**
     * Builds and returns the particle's 4-momentum as an array ordered as 
     * \f$(\mathbf p, E) = (p_x, p_y, p_z, E)\f$
     * @brief Returns the particle's 4-momentum
     * @return The particle's 4-momentum as a 4 components double array
     */
    /*inline double* P4() const {
      double out[4] = { fMomentum.P(0), fMomentum.P(1), fMomentum.P(2), E() };
      return out;
    }*/
    /*inline std::vector<double> P5() const {
      double* in = P4();
      std::vector<double> out(in, in+3);
      out.push_back(E());
      out.push_back(fMass);
      return out;
    }*/
    /**
     * @brief Sets the particle's energy
     * @param[in] E_ Energy, in GeV
     */
    inline void SetE(double E_=-1.) {
      if (E_<0.) E_ = sqrt(M2()+pow(P(),2));
      fMomentum.SetE(E_);
    }
    /**
     * @brief Gets the particle's energy in GeV
     */
    inline double E() const { return (fMomentum.E()<0.) ? std::sqrt(M2()+std::pow(P(),2)) : fMomentum.E(); };
    /**
     * @brief Gets the particle's squared energy in \f$\text{GeV}^\text{2}\f$
     */
    inline double E2() const { return std::pow(E(), 2); };
    void RotateThetaPhi(double theta_, double phi_);
    inline double Theta() const { return atan2(Pt(), Pz()); }
    inline double Phi() const { return atan2(Py(), Px()); }
    /**
     * @brief Is this particle a valid particle which can be used for kinematic computations ?
     */
    bool Valid();
    /**
     * Sets the "mother" of this particle (particle from which this particle is issued)
     * @brief Sets the mother particle (from which this particle arises)
     * @param[in] part_ A Particle object containing all the information on the mother particle
     */
    void SetMother(Particle* part_);
    /**
     * @brief Gets the unique identifier to the mother particle from which this particle arises
     * @return An integer representing the unique identifier to the mother of this particle in the event
     */
    inline ParticlesIds GetMothersIds() const { return fMothers; }
    /**
     * Adds a "daughter" to this particle (one of its decay product(s))
     * @brief Specify a decay product for this particle
     * @param[in] part_ The Particle object in which this particle will desintegrate or convert
     * @return A boolean stating if the particle has been added to the daughters list or if it was already present before
     */
    bool AddDaughter(Particle* part_);
    /**
     * @brief Gets the number of daughter particles arising from this one
     */
    inline unsigned int NumDaughters() { return fDaughters.size(); };
    /**
     * @brief Gets a vector containing all the daughters unique identifiers from this particle
     * @return An integer vector containing all the daughters' unique identifier in the event
     */
    std::vector<int> GetDaughters();
    void PDF2PDG();
    /**
     * Hadronises the particle with Pythia, and builds the shower (list of Particle objects) embedded in this object
     * @param[in] algo_ Algorithm in use to hadronise the particle
     * @brief Hadronises the particle using Pythia
     * @return A boolean stating whether or not the particle has been hadronised
     */
    bool Hadronise(std::string algo_);
    /**
     * @brief Is this particle a primary particle ?
     */
    inline bool Primary() const { return fIsPrimary; }
    /**
     * Codes 1-10 correspond to currently existing partons/particles, and larger codes contain partons/particles which no longer exist, or other kinds of event information
     * @brief Particle status
     */
    Status status;
  private:
    /**
     * @brief Mass in \f$\text{GeV}/c^2\f$
     */
    double fMass;
    /**
     * @brief List of mother particles
     */
    ParticlesIds fMothers;
    /**
     * @brief List of daughter particles
     */
    ParticlesIds fDaughters;
    /**
     * Unique identifier for a particle type. From @cite Beringer:1900zz :
     * _The Monte Carlo particle numbering scheme [...] is intended to facilitate interfacing between event generators, detector simulators, and analysis packages used in particle physics._
     * @brief Particle Data Group integer identifier
     */
    ParticleCode fPDGid;
    /**
     * @brief Is the particle a primary particle ?
     */
    bool fIsPrimary;
    Momentum fMomentum;
    double __tmp3[3], __tmp4[4];
};

inline bool compareParticle(Particle a, Particle b) { return a.id<b.id; }
inline bool compareParticlePtrs(Particle* a, Particle* b) { return a->id<b->id; }

typedef std::vector<Particle> Particles;
typedef std::vector<Particle*> ParticlesRef;
typedef std::multimap<int,Particle> ParticlesMap;

#endif
