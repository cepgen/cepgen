#ifndef _PARTICLE_H
#define _PARTICLE_H

#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <set>
#include <vector>
#include <algorithm>

#include "utils.h"
//#include "lheutils.h"

typedef std::set<int> ParticlesIds;

/**
 * Kinematic information for one particle
 * @brief Kinematics of one particle
 */
class Particle {
  public:
    enum ParticleCode {
      QUARK_D = 1,
      QUARK_U = 2,
      ELECTRON = 11,
      MUON = 13,
      TAU = 15,
      GLUON = 21,
      PHOTON = 22,
      PI_PLUS = 211,
      PI_0 = 111,
      RHO_770_0 = 113,
      OMEGA_782 = 223,
      J_PSI = 443,
      PHI_1680 = 100333,
      UPSILON_1S = 553,
      UPSILON_2S = 100553,
      UPSILON_3S = 200553,
      DIQUARK_UD0 = 2101,
      DIQUARK_UD1 = 2103,
      DIQUARK_UU1 = 2203,
      PROTON = 2212,
      NEUTRON = 2112,
      POMERON = 990,
      REGGEON = 110
    };
    /**
     * Gets the mass in GeV/c**2 of a particle given its PDG identifier
     * @brief Gets the mass of a particle
     * @param pdgId_ PDG ID of the particle
     * @return Mass of the particle in \f$\text{GeV}/c^2\f$
     */
    static double GetMassFromPDGId(Particle::ParticleCode pdgId_);
    /**
     * Gets the total decay width for one particle to be decayed
     * @param[in] pdgId_ PDG ID of the particle
     * @return Decay width in GeV
     */
    static double GetWidthFromPDGId(Particle::ParticleCode pdgId_);
    
    Particle();
    /**
     * @brief Object constructor (providing the role of the particle in the process, and its Particle Data Group identifier)
     */
    Particle(int role_,ParticleCode pdgId_=(ParticleCode)0);
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
    inline bool operator<(const Particle& rhs) { std::cout << this->id << "\t" << rhs.id << "\t" << (this->id<rhs.id) << std::endl; return this->id<rhs.id; };
    /**
     * @brief Comparison operator to enable the sorting of Particle objects' pointers in an event according to their reference's unique identifier
     */
    inline bool operator<(const Particle *rhs) { std::cout << this->id << "\t" << rhs->id << "\t" << (this->id<rhs->id) << std::endl; return this->id<rhs->id; };
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
    void LorentzBoost(double m_, double p_[4]);
    /**
     * Lorentz boost from ROOT
     */
    double* LorentzBoost(double p_[3]);
    /**
     * @brief Unique identifier of the particle (in a Event object context)
     */
    int id;
    /**
     * Unique identifier for a particle type. From @cite Beringer:1900zz :
     * _The Monte Carlo particle numbering scheme [...] is intended to facilitate interfacing between event generators, detector simulators, and analysis packages used in particle physics._
     * @brief Particle Data Group integer identifier
     */
    ParticleCode pdgId;
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
    /**
     * Particle's helicity
     * @fixme Float??
     */
    float helicity;
    /**
     * @brief Momentum along the \f$x\f$-axis in \f$\text{GeV}/c\f$
     */
    inline double Px() const { return this->_p4[0]; };
    /**
     * @brief Momentum along the \f$y\f$-axis in \f$\text{GeV}/c\f$
     */
    inline double Py() const { return this->_p4[1]; };
    /**
     * @brief Momentum along the \f$z\f$-axis in \f$\text{GeV}/c\f$
     */
    inline double Pz() const { return this->_p4[2]; };
    /**
     * Gets the particle's mass in \f$\text{GeV}/c^{2}\f$.
     * @brief Gets the particle's mass
     * @return The particle's mass
     */
    inline double M() const { return this->_m; };
    /**
     * Set the mass of the particle in \f$\text{GeV}/c^{2}\f$ according to a value given as an argument. This method ensures that the kinematics is properly set (the mass is set according to the energy and the momentum in priority)
     * @param m_ The mass in \f$\text{GeV}/c^{2}\f$ to set
     * @brief Set the particle's mass in \f$\text{GeV}/c^{2}\f$
     * @return A boolean stating whether or not the mass was correctly set
     */
    bool M(double m_);
    /**
     * @brief Gets the particle's squared mass
     */
    inline double M2() const { return std::pow(this->_m,2); };
    /**
     * Computes and returns \f$\eta\f$, the pseudo-rapidity of the particle
     * @brief Pseudo-rapidity
     * @return The pseudo-rapidity of the particle
     */
    inline double Eta() {
      return (this->Pt()!=0.)
	? log((this->P()+fabs(this->Pz()))/this->Pt())*(this->Pz()/fabs(this->Pz()))
	: 9999.*(this->Pz()/fabs(this->Pz()));
    };
    /**
     * Computes and returns \f$y\f$, the rapidity of the particle
     * @brief Rapidity
     * @return The rapidity of the particle
     */
    inline double Rapidity() { return (this->E()<0.) ? 999. : log((this->E()+this->Pz())/(this->E()-this->Pz()))/2.; };
    /**
     * Computes and returns \f$\phi\f$, the azimuthal angle of the particle in the transverse plane
     * @brief Azimuthal angle
     * @return The azimuthal angle of the particle
     */
    /*inline double Phi() const {
      return (this->Px()==0. && this->Py()==0. && this->Pz()==0.)
	? 0.
	: atan2(this->P(), this->Pz());
	};*/
    /**
     * @brief Sets the 3-momentum associated to the particle
     * @param[in] px_ Momentum along the \f$x\f$-axis, in \f$\text{GeV}/c\f$
     * @param[in] py_ Momentum along the \f$y\f$-axis, in \f$\text{GeV}/c\f$
     * @param[in] pz_ Momentum along the \f$z\f$-axis, in \f$\text{GeV}/c\f$
     * @return A boolean stating the validity of this particle (according to its 4-momentum norm)
     */
    inline bool P(double px_,double py_,double pz_) {
      double pp4[] = { px_, py_, pz_, -1. }; std::copy(pp4, pp4+4, this->_p4);
      this->M(-1.); this->E(-1.);
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
    inline bool P(double px_,double py_,double pz_,double E_) { double pp4[] = { px_, py_, pz_, E_ }; return this->P(pp4); };
    /**
     * @brief Sets the 4-momentum associated to the particle
     * @param[in] p_ 4-momentum
     * @return A boolean stating the validity of the particle's kinematics
     */
    inline bool P(double p_[4]) { 
      std::copy(p_, p_+4, this->_p4);
      if (p_[3]<0.) return this->P(p_[0], p_[1], p_[2]);
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
      if (c_>=0 and c_<4) return _p4[c_];
      else if (c_==4) return this->M();
      else return -999.;
    };
    /**
     * @brief Transverse momentum, in \f$\text{GeV}/c\f$
     */
    inline double Pt() const { return std::sqrt(std::pow(this->Px(),2)+std::pow(this->Py(),2)); };
    /**
     * @brief Norm of the 3-momentum, in \f$\text{GeV}/c\f$
     * @return The particle's 3-momentum norm as a double precision float
     */
    inline double P() const { return std::sqrt(std::pow(this->Pt(),2)+std::pow(this->Pz(),2)); };
    /**
     * Builds and returns the particle's 4-momentum as an array ordered as 
     * \f$(\mathbf p, E) = (p_x, p_y, p_z, E)\f$
     * @brief Returns the particle's 4-momentum
     * @return The particle's 4-momentum as a 4 components double array
     */
    inline double* P4() {
      this->_p4[3] = this->E();
      return this->_p4;
    };
    inline std::vector<double> P5() const {
      std::vector<double> out(this->_p4, this->_p4+3);
      out.push_back(this->E());
      out.push_back(this->_m);
      return out;
    }
    /**
     * @brief Sets the particle's energy
     * @param[in] E_ Energy, in GeV
     */
    inline void E(double E_) { this->_p4[3] = (E_<0.) ? std::sqrt(this->M2()+std::pow(this->P(),2)) : E_; };
    /**
     * @brief Gets the particle's energy in GeV
     */
    inline double E() const { return (this->_p4[3]<0.) ? std::sqrt(this->M2()+std::pow(this->P(),2)) : this->_p4[3]; };
    /**
     * @brief Gets the particle's squared energy in \f$\text{GeV}^\text{2}\f$
     */
    inline double E2() const { return std::pow(this->E(), 2); };
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
    inline ParticlesIds GetMothersIds() const { return this->_moth; }
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
    inline unsigned int NumDaughters() { return this->_daugh.size(); };
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
    inline bool Primary() const { return this->_isPrimary; }
    /**
     * Codes 1-10 correspond to currently existing partons/particles, and larger codes contain partons/particles which no longer exist, or other kinds of event information
     * @brief Particle status
     */
    int status;
  private:
    /**
     * @brief Mass in \f$\text{GeV}/c^2\f$
     */
    double _m;
    /**
     * @brief List of mother particles
     */
    ParticlesIds _moth;
    /**
     * @brief List of daughter particles
     */
    ParticlesIds _daugh;
    /**
     * @brief Is the particle a primary particle ?
     */
    bool _isPrimary;
    /**
     * List of components to characterise the particle's kinematics :
     * - 0-2: \f$\mathbf p = (p_x, p_y, p_z)\f$ (in \f$\text{GeV}/c\f$)
     * - 3: \f$E\f$ (in \f$\text{GeV}\f$)
     */
    double _p4[4];
    double __tmp3[3], __tmp4[4];
};

inline bool compareParticle(Particle a, Particle b) { return a.id<b.id; }
inline bool compareParticlePtrs(Particle* a, Particle* b) { return a->id<b->id; }

#endif
