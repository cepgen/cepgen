#ifndef _PARTICLE_H
#define _PARTICLE_H

#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <set>
#include <vector>

#include "utils.h"

/**
 * Kinematic information for one particle
 * @brief Kinematics of one particle
 */
class Particle {
  public:
    Particle();
    /**
     * @brief Object constructor (providing the role of the particle in the process, and its Particle Data Group identifier)
     */
    Particle(int role_,int pdgId_=0);
    ~Particle();
    /**
     * @brief Copies all the relevant quantities from one Particle object to another
     */
    Particle& operator=(const Particle&);
    /**
     * Returns a string containing all the particle's kinematics as expressed in the Les Houches format
     * @param revert_ Is the event symmetric ? If set to true, the third component of the momentum is reverted.
     * @return The LHE line associated to the particle, and containing the particle's history (mother/daughters), its kinematics, and its status
     */
    std::string GetLHEline(bool revert_=false);
    /**
     * Dumps into the standard output stream all the available information on this particle
     * @brief Dumps all the information on this particle
     */
    void Dump();
    /**
     * @brief Unique identifier of the particle (in a Event object context)
     */
    int id;
    /**
     * Unique identifier for a particle type. From @cite Beringer:1900zz :
     * _The Monte Carlo particle numbering scheme [...] is intended to facilitate interfacing between event generators, detector simulators, and analysis packages used in particle physics._
     * @brief Particle Data Group integer identifier
     */
    int pdgId;
    /**
     * @brief Role in the considered process
     */
    int role;
    /**
     * @brief Momentum along the \f$x\f$-axis in GeV/c
     */
    double px;
    /**
     * @brief Momentum along the \f$y\f$-axis in GeV/c
     */
    double py;
    /**
     * @brief Momentum along the \f$z\f$-axis in GeV/c
     */
    double pz;
    /**
     * Codes 1-10 correspond to currently existing partons/particles, and larger codes contain partons/particles which no longer exist, or other kinds of event information
     * @brief Particle status
     */
    int status;
    /**
     * Gets the particle's mass in GeV/c\f$^{2}\f$.
     * @brief Gets the particle's mass
     * @return The particle's mass
     */
    inline double M() { return this->m; };
    /**
     * @brief Set the particle's mass in GeV/c\f${}^\mathrm{2}\f$
     */
    bool M(double m_);
    /**
     * @brief Gets the particle's squared mass
     */
    inline double M2() { return std::pow(this->m,2); };
    /**
     * @brief Sets the 3-momentum associated to the particle
     * @param px_ Momentum along the \f$x\f$-axis, in GeV/c
     * @param py_ Momentum along the \f$y\f$-axis, in GeV/c
     * @param pz_ Momentum along the \f$z\f$-axis, in GeV/c
     * @return A boolean stating the validity of this particle (according to its 4-momentum norm)
     */
    inline bool P(double px_,double py_,double pz_) {
      this->px=px_; this->py=py_; this->pz=pz_;
      this->p3[0] = px_; this->p3[1] = py_; this->p3[2] = pz_;
      this->M(-1.);
      if (this->E()<0.) {
	if (this->M()>=0.) this->E(sqrt(pow(this->P(),2)-this->M2()));
	else return false;
      }
      else this->E(-1.);
      return true;
    };
    /**
     * Computes and returns \f$\eta\f$, the pseudo-rapidity of the particle
     * @brief Pseudo-rapidity
     * @return The pseudo-rapidity of the particle
     */
    inline double Eta() {
      return (this->Pt()!=0.)
	? log((this->P()+fabs(this->pz))/this->Pt())*(this->pz/fabs(this->pz))
	: 9999.*(this->pz/fabs(this->pz));
    }
    /**
     * Sets the 4-momentum associated to the particle, and computes its (invariant) mass.
     * @brief Sets the 4-momentum associated to the particle
     * @param px_ Momentum along the \f$x\f$-axis, in GeV/c
     * @param py_ Momentum along the \f$y\f$-axis, in GeV/c
     * @param pz_ Momentum along the \f$z\f$-axis, in GeV/c
     * @param E_ Energy, in GeV
     * @return A boolean stating the validity of the particle's kinematics
     */
    inline bool P(double px_,double py_,double pz_,double E_) {
      this->P(px_,py_,pz_);
      if (pow(E_,2)<pow(this->P(),2)) return false;
      if (this->M()<0.) this->M(sqrt(pow(this->E(),2)-pow(this->P(),2)));
      return true;
    };
    /**
     * @brief Sets the 4-momentum associated to the particle
     * @param p_ 3-momentum
     * @param E_ Energy, in GeV
     * @return A boolean stating the validity of the particle's kinematics
     */
    bool P(double p_[3],double E_);
    /**
     * @brief Sets the 4-momentum associated to the particle
     * @param p_ 4-momentum
     * @return A boolean stating the validity of the particle's kinematics
     */
    inline bool P(double p_[4]) { 
      double p3[3];
      std::copy(p_, p_+3, p3);
      return this->P(p3, p_[3]);
    };
    /**
     * @brief Returns the particle's 3-momentum
     * @return The particle's 3-momentum as a 3 components double array
     */
    inline double* P3() { return this->p3; };
    /**
     * Builds and returns the particle's 4-momentum as an array ordered as 
     * \f$(\mathbf p, E) = (p_x, p_y, p_z, E)\f$
     * @brief Returns the particle's 4-momentum
     * @return The particle's 4-momentum as a 4 components double array
     */
    inline double* P4() {
      std::copy(this->p3, this->p3+3, this->p4);
      this->p4[3] = this->E();
      return this->p4;
    };
    /**
     * @brief Norm of the 3-momentum, in GeV/c
     * @return The particle's 3-momentum norm as a double precision float
     */
    inline double P() { return std::sqrt(std::pow(this->Pt(),2)+std::pow(this->pz,2)); };
    /**
     * @brief Transverse momentum, in GeV/c
     */
    inline double Pt() { return std::sqrt(std::pow(this->px,2)+std::pow(this->py,2)); };
    /**
     * @brief Sets the particle's energy
     * @param E_ Energy, in GeV
     */
    inline void E(double E_) { this->e=E_; };
    /**
     * @brief Gets the particle's energy
     */
    inline double E() { return this->e; };
    /**
     * @brief Is this particle a valid particle which can be used for kinematic computations ?
     */
    bool Valid();
    /**
     * Sets the "mother" of this particle (particle from which this particle is issued)
     * @brief Sets the mother particle (from which this particle arises)
     * @param part_ A Particle object containing all the information on the mother particle
     */
    void SetMother(Particle* part_);
    /**
     * @brief Gets the mother particle from which this particle arises
     */
    Particle* GetMother();
    /**
     * Adds a "daughter" to this particle (one of its decay product(s))
     * @brief Specify a decay product for this particle
     * @param part_ The Particle object in which this particle will desintegrate or convert
     * @return A boolean stating if the particle has been added to the daughters list or if it was already present before
     */
    bool AddDaughter(Particle* part_);
    /**
     * @brief Gets the number of daughter particles arising from this one
     */
    inline unsigned int NumDaughters() { return this->_daugh.size(); };
    /**
     * @brief Gets a vector containing all the daughters from this particle
     * @return A Particle objects vector containing all the daughters' kinematic information
     */
    std::vector<Particle*> GetDaughters();
    void PDF2PDG();
    /**
     * Hadronises the particle with Pythia, and builds the shower (list of Particle objects) embedded in this object
     * @param algo_ Algorithm in use to hadronise the particle
     * @brief Hadronises the particle using Pythia
     */
    bool Hadronise(std::string algo_);
  private:
    /**
     * @brief Energy, in GeV
     */
    double e;
    /**
     * @brief Mass in GeV/c\f$^{2}\f$
     */
    double m;
    /**
     * @brief Mother particle
     */
    Particle *_mother;
    /**
     * @brief List of daughter particles
     */
    std::set<Particle*> _daugh;
    /**
     * @brief Is the particle a primary particle ?
     */
    bool _isPrimary;
    double p3[3];
    double p4[4];
};

#endif
