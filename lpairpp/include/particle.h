#ifndef _PARTICLE_H
#define _PARTICLE_H

#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <vector>

#include "utils.h"

/**
 * Kinematic information for one particle
 * @brief Kinematics of one particle
 */
class Particle {
  public:
    Particle();
    Particle(int role_,int pdgId_=0);
    ~Particle();
    /**
     * Returns a string containing all the particle's kinematics as expressed in
     * the Les Houches format
     * @param revert_ Is the event symmetric ? If set to true, the third component
     * of the momentum is reverted.
     * @return The LHE line associated to the particle, and containing the particle's
     * history (mother/daughters), its kinematics, and its status
     */
    std::string GetLHEline(bool revert_=false);
    /**
     * Dumps into the standard output stream all the available information on this
     * particle
     * @brief Dumps all the information on this particle
     */
    void Dump();
    /** @brief Particle Data Group integer identifier */
    int pdgId;
    /** @brief Role in the considered process */
    int role;
    /** @brief Momentum along the \f$x\f$-axis in GeV/c */
    double px;
    /** @brief Momentum along the \f$y\f$-axis in GeV/c */
    double py;
    /** @brief Momentum along the \f$z\f$-axis in GeV/c */
    double pz;
    /** @brief Transverse momentum, in GeV/c */
    double pt;
    /** @brief Norm of the 3-momentum, in GeV/c */
    double p;
    /** @brief Pseudo-rapidity */
    double eta;
    /**
     * Codes 1-10 correspond to currently existing partons/particles, and larger
     * codes contain partons/particles which no longer exist, or other kinds of
     * event information
     * @brief Particle status
     */
    int status;
    /**
     * Gets the particle's mass in GeV/c\f$^{2}\f$.
     * @brief Gets the particle's mass
     * @return The particle's mass
     */
    inline double M() { return this->m; };
    /** @brief Set the particle's mass in GeV/c\f${}^\mathrm{2}\f$ */
    bool M(double m_);
    /**
     * @brief Sets the 3-momentum associated to the particle
     * @param px_ Momentum along the \f$x\f$-axis, in GeV/c
     * @param py_ Momentum along the \f$y\f$-axis, in GeV/c
     * @param pz_ Momentum along the \f$z\f$-axis, in GeV/c
     */
    inline bool P(double px_,double py_,double pz_) {
      this->px=px_; this->py=py_; this->pz=pz_;
      this->pt=sqrt(pow(px_,2)+pow(py_,2));
      this->p=sqrt(pow(this->pt,2)+pow(pz_,2));
      this->M(-1.);
      if (this->E()<0.) {
	if (this->M()>=0.) this->E(sqrt(pow(this->p,2)-pow(this->M(),2)));
	else return false;
      }
      else this->E(-1.);
      this->eta=(this->pt!=0.)
	? log((this->p+fabs(this->pz))/this->pt)*(this->pz/fabs(this->pz))
	: 9999.*(this->pz/fabs(this->pz));
      //this->eta=0.5*log((this->p+this->pz)/(this->p-this->pz));
      return true;
    };
    /**
     * Sets the 4-momentum associated to the particle, and 
     * computes its (invariant) mass.
     * @brief Sets the 4-momentum associated to the particle
     * @param px_ Momentum along the \f$x\f$-axis, in GeV/c
     * @param py_ Momentum along the \f$y\f$-axis, in GeV/c
     * @param pz_ Momentum along the \f$z\f$-axis, in GeV/c
     * @param E_ Energy, in GeV
     */
    inline bool P(double px_,double py_,double pz_,double E_) {
      this->P(px_,py_,pz_);
      if (pow(E_,2)<pow(this->p,2)) {
        return false;
      }
      if (this->M()<0.) {
	this->M(sqrt(pow(this->E(),2)-pow(this->p,2)));
      }
      return true;
    };
    /**
     * @brief Sets the 4-momentum associated to the particle
     * @param p_ 3-momentum
     * @param E_ Energy, in GeV
     */
    bool P(double p_[3],double E_);
    /**
     * @brief Sets the particle's energy
     * @param E_ Energy, in GeV
     */
    inline void E(double E_) { this->e=E_; };
    /** @brief Gets the particle's energy */
    inline double E() { return this->e; };
    /** @brief Gets the particle's squared mass */
    inline double M2() { return std::pow(this->m,2); };
    /** @brief Is this particle a valid particle which can be used for kinematic computations ? */
    bool Valid();
    /**
     * @brief Sets the mother particle (from which this particle arises)
     * @param part_ A Particle object containing all the information on the mother particle
     */
    inline void SetMother(Particle* part_) { this->_mother=part_; this->_isPrimary=false; };
    /**
     * @brief Gets the mother particle from which this particle arises
     */
    Particle* GetMother();
    /**
     * @brief Specify a decay product for this particle
     * @param part_ The Particle object in which this particle will desintegrate or convert
     */
    void AddDaughter(Particle* part_);
    /** @brief Gets the number of daughter particles arising from this one */
    inline unsigned int NumDaughters() { return this->_daugh->size(); };
    /**
     * @brief Gets a daughter from this particle, labelled by its identifier in this particle's daughters list
     * @return A Particle object containing all the kinematic information related to this daughter particle
     */
    Particle* GetDaughter(const unsigned int num_=0);
  private:
    /** @brief Energy, in GeV */
    double e;
    /** @brief Mass in GeV/c\f$^{2}\f$ */
    double m;
    Particle *_mother;
    std::vector<Particle> *_daugh;
    bool _isPrimary;
};

#endif
