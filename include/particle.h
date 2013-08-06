#ifndef _PARTICLE_H
#define _PARTICLE_H

#include <iostream>
#include <sstream>
#include <cmath>
#include <string>

/**
 * Kinematic information for one particle
 * @brief Kinematics of one particle
 */
class Particle {
  public:
    Particle();
    ~Particle();
    /**
     * Returns a string containing all the particle's kinematics as expressed in
     * the Les Houches format
     * @param revert_ Is the event symmetric ? If set to true, the third component
     * of the momentum is reverted.
     * @return The LHE line
     */
    std::string GetLHEline(bool revert_=false);
    /** @brief Particle Data Group integer identifier */
    int pdgId;
    /** @brief Role in the considered process */
    int role;
    /** @brief Energy in GeV */
    double e;
    /** @brief Mass in GeV/c\f$^{2}\f$ */
    double m;
    /** @brief Momentum along the \f$x\f$-axis in GeV/c */
    double px;
    /** @brief Momentum along the \f$y\f$-axis in GeV/c */
    double py;
    /** @brief Momentum along the \f$z\f$-axis in GeV/c */
    double pz;
    /** @brief Transverse momentum */
    double pt;
    /** @brief Momentum */
    double p;
    /**
     * @brief Sets the 3-momentum associated to the particle
     * @param px_ Momentum along the \f$x\f$-axis
     * @param py_ Momentum along the \f$y\f$-axis
     * @param pz_ Momentum along the \f$z\f$-axis
     */
    inline void SetP(double px_,double py_,double pz_) {
      this->px=px_; this->py=py_; this->pz=pz_; this->pt=sqrt(pow(px_,2)+pow(py_,2)); this->p=sqrt(pow(this->pt,2)+pow(pz_,2));
    };
    /**
     * @brief Sets the 4-momentum associated to the particle
     * @param px_ Momentum along the \f$x\f$-axis
     * @param py_ Momentum along the \f$y\f$-axis
     * @param pz_ Momentum along the \f$z\f$-axis
     * @param E_ Energy
     */
    inline void SetP(double px_,double py_,double pz_,double E_) {
      this->SetP(px_,py_,pz_); this->SetE(E_); this->m=(sqrt(pow(E_,2)-pow(this->p,2)));
    };
    inline void SetE(double E_) {
      this->e=E_;
    };
};

#endif
