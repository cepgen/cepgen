#ifndef _GAMGAM_H
#define _GAMGAM_H

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <map>

#include "event.h"
#include "utils.h"

/**
 * @brief List of kinematic cuts to apply on the central and outgoing phase space.
 */
class GamGamKinematics {
  public:
    GamGamKinematics();
    ~GamGamKinematics();
    void Dump();
    /**
     * Type of kinematics to consider for the process. Can either be :
     *  * 0 for the electron-electron elastic case
     *  * 1 for the proton-proton elastic case
     *  * 2 for the proton-proton single-dissociative (or inelastic) case
     *  * 3 for the proton-proton double-dissociative case
     * @brief Type of kinematics to consider for the phase space
     */
    int kinematics;
    /** @brief Sets of cuts to apply on the final phase space */
    int mode;
    /** @brief Minimal transverse momentum of the single outgoing leptons */
    double ptmin;
    /** @brief Maximal transverse momentum of the single outgoing leptons */
    double ptmax;
    /** @brief Minimal energy of the central two-photons system */
    double emin;
    /** @brief Maximal energy of the central two-photons system */
    double emax;
    /** @brief Minimal polar (\f$\theta_\mathrm{min}\f$) angle of the outgoing leptons, expressed in degrees */
    double thetamin;
    /** @brief Maximal polar (\f$\theta_\mathrm{max}\f$) angle of the outgoing leptons, expressed in degrees */
    double thetamax;
    double mxmin;
    double mxmax;
    /** @brief The minimal value of \f$Q^2\f$ */
    double q2min;
    /** @brief The maximal value of \f$Q^2\f$ */
    double q2max;
    /** @brief The minimal \f$s\f$ on which the cross section is integrated */
    double wmin;
    /** @brief The maximal \f$s\f$ on which the cross section is integrated.
     *  If negative, the maximal energy available to the system (hence,
     * \f$s=(\sqrt{s})^{2}\f$) is provided.*/
    double wmax;
};

/**
 * Full class of methods and objects to compute the full analytic matrix element
 * @cite Vermaseren1983347 for the \f$\gamma\gamma\to\ell^{+}\ell^{-}\f$ process
 * according to a set of kinematic constraints provided for the incoming and
 * outgoing particles (the GamGamKinematics object).
 * @brief Computes the matrix element for a \f$\gamma\gamma\to\ell^{+}\ell^{-}\f$
 *  process
 */
class GamGam {
 public:
  /**
   * Sets the mandatory parameters used in the methods computing the kinematics
   * and the cross-section of this phase space point.
   * @brief Class constructor
   * @param ndim_ The number of dimensions of the point in the phase space
   * @param nOpt_ Optimisation???
   * @param x_[] The (@a ndim_)-dimensional point in the phase space on which
   *  the kinematics and the cross-section are computed
   * @todo Figure out how this @a nOpt_ parameter is affecting the final
   *  cross-section computation and events generation
   */
  GamGam(const unsigned int ndim_,int nOpt_,double x_[]);
  ~GamGam();
  /**
   * Specifies the incoming particles' kinematics as well as their properties
   *  (role in the process and code according to the PDG convention)
   * @brief Sets the momentum and PDG id for the incoming particles
   * @param part_ Role of the particle in the process
   * @param momentum_[] 3-momentum of the particle
   * @param pdgId_ Particle ID according to the PDG convention
   * @return True if the kinematics was correctly set for the given particle role
   */
  bool SetIncomingKinematics(int part_,double momentum_[3],int pdgId_);
  /**
   * Specifies the incoming particles' kinematics as well as their properties
   * using two Particle objects.
   * @brief Sets the momentum and PDG id for the incoming particles
   * @param ip1_ Information on the first incoming particle
   * @param ip2_ Information on the second incoming particle
   */
  bool SetIncomingKinematics(Particle ip1_,Particle ip2_);
  /**
   * @brief Sets the PDG id for the outgoing particles
   * @param part_ Role of the particle in the process
   * @param pdgId_ Particle ID according to the PDG convention
   */
  bool SetOutgoingParticles(int part_,int pdgId_);
  /*void SetCutsMode(int mod_) {this->_modcut = mod_;};
  //void SetCutsMode(int);
  void SetMinimumPt(double ptmin_) {this->_ptcut = ptmin_;};
  void SetMinimumE(double emin_) {this->_ecut = emin_;};
  void SetThetaRange(double,double);*/
  /**
   * @brief Sets the list of kinematic cuts to apply on the outgoing particles'
   *  final state
   * @param cuts_ The Cuts object containing the kinematic parameters
   */
  void SetCuts(GamGamKinematics cuts_);
  /**
   * @brief Get a particle given its role in the process
   * @param role_ An integer denoting the particle's role in the selected
   *  production process
   */
  Particle* GetParticle(int role_);
  void SetParticle(int,Particle*);
  /**
   * Is the system's kinematics well defined and compatible with the process ?
   * This check is mandatory to perform the (@a _ndim)-dimensional point's
   * cross-section computation.
   * @brief Is the system's kinematics well defined?
   * @return A boolean stating if the input kinematics and the final
   *  states are well defined
   */
  bool IsKinematicsDefined() { return setkin; }
  /**
   * Computes the centre of mass energy for the system,
   *  according to the incoming particles' kinematics
   * @brief Computes \f$\sqrt{s}\f$ for the system
   */
  void ComputeSqS();
  /**
   * Computes the cross-section for the \f$\gamma\gamma\to\ell^{+}\ell^{-}\f$
   *  process with the given kinematics
   * @brief Computes the process' cross section
   * @return \f$\frac{\textrm d\sigma}{\mathrm d\mathbf x}(\gamma\gamma\to\ell^{+}\ell^{-})\f$,
   * the differential cross-section for the given point in the phase space.
   */
  double ComputeXsec(int nm_=1);
  void FillKinematics(bool symmetrise_=false);
  void StoreEvent(std::ofstream*,double);
  double GetT1() { return this->_t1; };
  void GetT1extrema(double& t1min_, double& t1max_) { t1min_=this->_t1min; t1max_=this->_t1max; };
  double GetT2() { return this->_t2; };
  void GetT2extrema(double& t2min_, double& t2max_) { t2min_=this->_t2min; t2max_=this->_t2max; };
  inline Event* GetEvent() { return this->_ev; };
 private:
  /**
   * Calculates energies and momenta of the 1st, 2nd (resp. the "proton-like"
   * and the "electron-like" incoming particles), 3rd (the "proton-like" outgoing
   * particle), 4th (the two-photons central system) and 5th (the "electron-like"
   * outgoing particle) particles in the overall centre of mass frame.
   * @brief Energies/momenta computation for the various particles, in the CM
   *  system
   */
  bool Orient();
  /**
   * Contains the expression of the matrix element squared for the process under
   * considerations. It returns the value of the convolution of the form factor
   * or structure functions with the central two-photons matrix element squared.
   * @brief Computes the matrix element squared for the requested process
   * @return The full matrix element for the two-photon production of a pair of
   *  spin\f$-\frac{1}{2}-\f$point particles
   */
  double PeriPP(int,int);
  /**
   * Describes the kinematics of the process \f$p_1+p_2\to p_3+p_4+p_5\f$ in
   * terms of Lorentz-invariant variables. These variables (along with others)
   * will then be feeded into the PeriPP method (thus are essential for the
   * evaluation of the full matrix element).
   */
  bool Pickin();
  /** @brief Number of dimensions on which the integration has to be performed  */
  unsigned int _ndim;
  /**
   * @brief Array of @a _ndim components representing the point on which the
   *  weight in the cross-section is computed
   */
  double *_x;
  int _nOpt;
  // COMMON/PICKZZ/
  /** @brief \f$\mathbf p_1\f$, 3-momentum of the first proton-like incoming particle */
  double _p3_p1[3];
  /** @brief \f$\left|\mathbf p_1\right|\f$, 3-momentum norm of the first proton-like incoming particle */
  double _pp1;
  /** @brief \f$E_1\f$, energy of the first proton-like incoming particle */
  double _ep1;
  /** @brief \f$m_1\f$, mass of the first proton-like incoming particle */
  double _mp1;
  /** @brief \f$m_1^2\f$, squared mass of the first proton-like incoming particle */
  double _w1;
  /** @brief PDG identifier of the first proton-like incoming particle */
  int _pdg1;
  /** @brief \f$\mathbf p_2\f$, 3-momentum of the second incoming particle */
  double _p3_p2[3];
  /** @brief \f$\left|\mathbf p_2\right|\f$, 3-momentum norm of the second proton-like incoming particle */
  double _pp2;
  /** @brief \f$E_2\f$, energy of the second proton-like incoming particle */
  double _ep2;
  /** @brief \f$m_2\f$, mass of the second proton-like incoming particle */
  double _mp2;
  /** @brief \f$m_2^2\f$, squared mass of the second proton-like incoming particle */
  double _w2;
  /** @brief PDG identifier of the second proton-like incoming particle */
  int _pdg2;
  /** @brief \f$\mathbf p_3\f$, 3-momentum of the first proton-like outgoing particle */
  double _p3_p3[3];
  /** @brief \f$\left|\mathbf p_3\right|\f$, 3-momentum norm of the first proton-like outgoing particle */
  double _pp3;
  /** @brief \f$E_3\f$, energy of the first proton-like outgoing particle */
  double _ep3;
  /** @brief \f$m_3\f$, mass of the first proton-like outgoing particle */
  double _mp3;
  /** @brief \f$m_3^2\f$, squared mass of the first proton-like outgoing particle */
  double _w3;
  /** @brief PDG identifier of the first proton-like outgoing particle */
  int _pdg3;
  /** @brief \f$\mathbf p_4\f$, 3-momentum of the two-photon central system */
  double _p3_c4[3];
  /** @brief \f$\left|\mathbf p_4\right|\f$, 3-momentum norm of the two-photon central system */
  double _pc4;
  /** @brief \f$E_4\f$, energy of the two-photon central system */
  double _ec4;
  /** @brief \f$m_4\f$, mass of the two-photon central system */
  double _mc4;
  /** @brief \f$m_4^2\f$, squared mass of the two-photon central system */
  double _w4;
  /** @brief \f$\mathbf p_5\f$, 3-momentum of the second proton-like outgoing particle */
  double _p3_p5[3];
  /** @brief \f$\left|\mathbf p_5\right|\f$, 3-momentum norm of the second proton-like outgoing particle */
  double _pp5;
  /** @brief \f$E_5\f$, energy of the second proton-like outgoing particle */
  double _ep5;
  /** @brief \f$m_5\f$, mass of the second proton-like outgoing particle */
  double _mp5;
  /** @brief \f$m_5^2\f$, squared mass of the second proton-like outgoing particle */
  double _w5;
  /** @brief PDG identifier of the second proton-like outgoing particle */
  int _pdg5;
  /** @brief \f$\mathbf p_6\f$, 3-momentum of the first outgoing lepton */
  double _p3_l6[3];
  /** @brief \f$\left|\mathbf p_6\right|\f$, 3-momentum norm of the first outgoing lepton */
  double _pl6;
  /** @brief \f$E_6\f$, energy of the first outgoing lepton */
  double _el6;
  /** @brief \f$m_6\f$, mass of the first outgoing lepton */
  double _ml6;
  /** @brief \f$m_6^2\f$, squared mass of the first outgoing lepton */
  double _w6;
  /** @brief \f$p_{T,6}\f$, transverse momentum of the first outgoing lepton */
  double _pt_l6;
  /** @brief \f$E_6^\mathrm{lab}\f$, energy of the first outgoing lepton, computed in the lab frame */
  double _e6lab;
  /** @brief PDG identifier of the first outgoing lepton */
  int _pdg6;
  /** @brief \f$\mathbf p_7\f$, 3-momentum of the second outgoing lepton */
  double _p3_l7[3];
  /** @brief \f$\left|\mathbf p_7\right|\f$, 3-momentum norm of the second outgoing lepton */
  double _pl7;
  /** @brief \f$E_7\f$, energy of the second outgoing lepton */
  double _el7;
  /** @brief \f$m_7\f$, mass of the second outgoing lepton */
  double _ml7;
  /** @brief \f$m_7^2\f$, squared mass of the second outgoing lepton */
  double _w7;
  /** @brief \f$p_{T,7}\f$, transverse momentum of the second outgoing lepton */
  double _pt_l7;
  /** @brief \f$E_7^\mathrm{lab}\f$, energy of the second outgoing lepton, computed in the lab frame */
  double _e7lab;
  /** @brief PDG identifier of the second outgoing lepton */
  int _pdg7;
  /** @brief Energy of the first central photon of momentum \f$t_1\f$ */
  double _eg1;
  /** @brief 3-momentum of the second central photon of momentum \f$t_1\f$ */
  double _p3_g1[3];
  /** @brief Energy of the second central photon of momentum \f$t_2\f$ */
  double _eg2;
  /** @brief 3-momentum of the second central photon of momentum \f$t_2\f$ */
  double _p3_g2[3];
  // CM energy of the incoming particles' system
  /** @brief \f$s\f$, squared centre of mass energy of the incoming particles' system */
  double _s;
  /** @brief \f$\sqrt{s}\f$, centre of mass energy of the incoming particles' system */
  double _sqs;
  /** @brief Total energy provided by the two incoming proton-like particles */
  double _etot;
  /** @brief Total momentum provided by the two incoming proton-like particles (along the \f$z\f$-axis) */
  double _ptot;
  /** @brief Minimal \f$Q^2\f$ exchange */
  double _q2min;
  /** @brief Maximal \f$Q^2\f$ exchange */
  double _q2max;
  double _qp2min, _qp2max;
  // COMMON /ACCURA/
  double _acc3;
  double _acc4;
  // COMMON /ANGU/
  /** @brief \f$\cos\theta_6^\mathrm{CM}\f$, production angle of the first outgoing lepton, computed in the centre of mass system. */
  double _ctcm6;
  /** @brief \f$\sin\theta_6^\mathrm{CM}\f$, production angle of the first outgoing lepton, computed in the centre of mass system. */
  double _stcm6;
  // COMMON /CIVITA/
  double _epsi;
  double _g5, _g6;
  double _a5, _a6;
  double _bb;
  // COMMON /DOTP/
  double _p12, _p13, _p14, _p15, _p23, _p24, _p25, _p34, _p35, _p45;
  double _p1k2, _p2k1;
  // COMMON /DOTPS/
  double _d1dq, _d1dq2, _q1dq, _q1dq2;
  // COMMON /EXTRA/
  double _s1, _s2;
  double _t1, _t1min, _t1max;
  double _t2, _t2min, _t2max;
  // COMMON /LEVI/
  double _gram;
  double _dd1, _dd2, _dd3, _dd5;
  double _delta;
  double _g4;
  double _sa1, _sa2;
  // COMMON /LTCOM/
  /**
   * @brief \f$\gamma\f$ factor of the centre of mass system, used in the
   * computation of the inverse boost for the outgoing leptons
   */
  double _gamma;
  /**
   * @brief \f$\beta\gamma\f$ factor of the centre of mass system, used in the
   * computation of the inverse boost for the outgoing leptons
   */
  double _betgam;
  // COMMON /LEVI/
  /**
   * @brief \f$\delta_1=m_3^2-m_1^2\f$ as defined in Vermaseren's paper
   * @cite Vermaseren1983347 for the full definition of this quantity
   */
  double _w31;
  /**
   * @brief \f$\delta_4=m_5^2-m_2^2\f$ as defined in Vermaseren's paper
   * @cite Vermaseren1983347 for the full definition of this quantity
   */
  double _w52;
  /**
   * @brief \f$\delta_5=m_4^2-t_1\f$ as defined in Vermaseren's paper
   * @cite Vermaseren1983347 for the full definition of this quantity
   */
  double _dd4;
  /**
   * @brief \f$\delta_2=m_1^2-m_2^2\f$ as defined in Vermaseren's paper
   * @cite Vermaseren1983347 for the full definition of this quantity
   */
  double _w12;
  /**
   * @brief \f$\delta_6=m_4^2-m_5^2\f$ as defined in Vermaseren's paper
   * @cite Vermaseren1983347 for the full definition of this quantity
   */
  double _tau;
  // COMMON /PICKZZ/
  double _sl1;
  // COMMON /QVEC/   // 0 = E, 1-3 = p
  double _qve[4];
  // COMMON /VARIAB/
  double _p;
  /** @brief \f$\cos\theta_3\f$ of the first outgoing proton-like particle */
  double _ct3;
  /** @brief \f$\sin\theta_3\f$ of the first outgoing proton-like particle */
  double _st3;
  /** @brief \f$\cos\theta_4\f$ of the two-photons centre of mass system */
  double _ct4;
  /** @brief \f$\sin\theta_4\f$ of the two-photons centre of mass system */
  double _st4;
  /** @brief \f$\cos\theta_5\f$ of the second outgoing proton-like particle */
  double _ct5;
  /** @brief \f$\sin\theta_5\f$ of the second outgoing proton-like particle */
  double _st5;
  /** @brief \f$\cos\phi_3\f$ of the first outgoing proton-like particle */
  double _cp3;
  /** @brief \f$\sin\phi_3\f$ of the first outgoing proton-like particle */
  double _sp3;
  /** @brief \f$\cos\phi_5\f$ of the second outgoing proton-like particle */
  double _cp5;
  /** @brief \f$\sin\phi_5\f$ of the second outgoing proton-like particle */
  double _sp5;
  // COMMON /VARIAC/
  double _al3, _al4;
  double _be4, _be5;
  double _de3, _de5;
  double _p_p3, _p_p4, _p_p5;
  // COMMON /VARIAD/
  /** @brief \f$\cos\theta_6\f$ of the first outgoing lepton */
  double _ct6;
  /** @brief \f$\sin\theta_6\f$ of the first outgoing lepton */
  double _st6;
  /** @brief \f$\cos\theta_7\f$ of the second outgoing lepton */
  double _ct7;
  /** @brief \f$\sin\theta_7\f$ of the second outgoing lepton */
  double _st7;
  /** @brief \f$\cos\phi_6\f$ of the first outgoing lepton */
  double _cp6;
  /** @brief \f$\sin\phi_6\f$ of the first outgoing lepton */
  double _sp6;
  /** @brief \f$\cos\phi_7\f$ of the second outgoing lepton */
  double _cp7;
  /** @brief \f$\sin\phi_7\f$ of the second outgoing lepton */
  double _sp7;
  double _dj;
  /** @brief Is the first incoming proton-like particle's kinematic set ? */
  bool setp1;
  /** @brief Is the second incoming proton-like particle's kinematic set ? */
  bool setp2;
  /** @brief Is the first outgoing proton-like particle's kinematic set ? */
  bool setp3;
  /** @brief Is the second outgoing proton-like particle's kinematic set ? */
  bool setp5;
  /** @brief Is the outgoing leptons' state set ? */
  bool setll;
  /** @brief Is the incoming particles' kinematic set ? */
  bool setin;
  /** @brief Is the outgoing particles' kinematic set ? */
  bool setout;
  /** @brief Is the full event's kinematic set ? */
  bool setkin;

  // Cuts
  /** @brief Set of cuts to apply on the final phase space */
  GamGamKinematics _cuts;
  double _cotth1, _cotth2;
  // Particles
  Event *_ev;
  //std::map<int,Particle> *_part;
};

#endif

