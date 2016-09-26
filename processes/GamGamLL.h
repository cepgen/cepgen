#ifndef GamGamLL_h
#define GamGamLL_h

#define PROCESS_DESCRIPTION "Two-photon production of lepton pairs"

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <map>

#include "core/GenericProcess.h"

/**
 * Full class of methods and objects to compute the full analytic matrix element
 * \cite Vermaseren1983347 for the \f$\gamma\gamma\to\ell^{+}\ell^{-}\f$ process
 * according to a set of kinematic constraints provided for the incoming and
 * outgoing particles (the Kinematics object).
 * The particle roles in this process are defined as following : \n
 * \image latex lpair_kinematics.pdf Detailed particle roles in the two-photon process as defined by the @a GamGamLL object. The incoming protons/electrons are denoted by a role 1, and 2, as the outgoing protons/protons remnants/ electrons carry the indices 3 and 5. The two outgoing leptons have the roles 6 and 7, while the lepton/antilepton distinction is done randomly (thus, the arrow convention is irrelevant here).
 *
 * The \a f function created by this Process child has its \a _ndim -dimensional
 * coordinates mapped as :
 * - 0 = \f$t_1\f$, first incoming photon's virtuality
 * - 1 = \f$t_2\f$, second incoming photon's virtuality
 * - 2 = \f$s_2\f$ mapping
 * - 3 = yy4 = \f$\cos\left(\pi x_3\right)\f$ definition
 * - 4 = \f$w_4\f$, the two-photon system's invariant mass
 * - 5 = xx6 = \f$\frac{1}{2}\left(1-\cos\theta^\text{CM}_6\right)\f$ definition (3D rotation of the first outgoing lepton with respect to the two-photon centre-of-mass system). If the @a nm_ optimisation flag is set this angle coefficient value becomes 
 *   \f[\frac{1}{2}\left(\frac{a_\text{map}}{b_\text{map}}\frac{\beta-1}{\beta+1}+1\right)\f]
 *   with \f$a_\text{map}=\frac{1}{2}\left(w_4-t_1-t_2\right)\f$, \f$b_\text{map}=\frac{1}{2}\sqrt{\left(\left(w_4-t_1-t_2\right)^2-4t_1t_2\right)\left(1-4\frac{w_6}{w_4}\right)}\f$, and \f$\beta=\left(\frac{a_\text{map}+b_\text{map}}{a_\text{map}-b_\text{map}}\right)^{2x_5-1}\f$
 *   and the \a fJacobian element is scaled by a factor \f$\frac{1}{2}\frac{\left(a_\text{map}^2-b_\text{map}^2\cos^2\theta^\text{CM}_6\right)}{a_\text{map}b_\text{map}}\log\left(\frac{a_\text{map}+b_\text{map}}{a_\text{map}-b_\text{map}}\right)\f$
 * - 6 = _phicm6_, or \f$\phi_6^\text{CM}\f$ the rotation angle of the dilepton system in the centre-of-mass
 *   system
 * - 7 = \f$x_q\f$, \f$w_X\f$ mappings, as used in the single- and double-dissociative
 *   cases only
 * \brief Compute the matrix element for a CE \f$\gamma\gamma\to\ell^{+}\ell^{-}\f$
 *  process
 */
class GamGamLL : public GenericProcess
{
 public:
  /**
   * Set the mandatory parameters used in the methods computing the kinematics
   * and the cross-section for the \f$\gamma\gamma\rightarrow\ell^{+}\ell^{-}\f$
   * process
   * \brief Class constructor
   * \param[in] nOpt_ Optimisation (legacy from LPAIR)
   */
  GamGamLL(int nOpt_=0);
  
  //~GamGamLL();
  void AddEventContent();
  void BeforeComputeWeight();
  /**
   * Compute the cross-section for the \f$\gamma\gamma\to\ell^{+}\ell^{-}\f$
   * process with the given kinematics
   * \brief Computes the process' weight for the given point
   * \return \f$\mathrm d\sigma(\mathbf x)(\gamma\gamma\to\ell^{+}\ell^{-})\f$,
   * the differential cross-section for the given point in the phase space.
   */
  double ComputeWeight();
  int GetNdim(Kinematics::ProcessMode) const;
  void FillKinematics(bool);
  /// Compute the ougoing proton remnant mass
  /// \param[in] x_ A random number (between 0 and 1)
  /// \param[in] outmass_ The maximal outgoing particles' invariant mass
  /// \param[in] lepmass_ The outgoing leptons' mass
  /// \param[out] dw_ The size of the integration bin
  /// \return Mass of the outgoing proton remnant
  double ComputeOutgoingPrimaryParticlesMasses(double x_, double outmass_, double lepmass_, double* dw_);
  /// Return the value of the first inner photon's virtuality
  /// \return \f$t_1\f$, the first photon virtuality
  inline double GetT1() const { return this->fT1; };
  /// Return the two limit values for the first photon's virtuality
  /// \param[out] t1min_ The minimal value for \f$t_1\f$
  /// \param[out] t1max_ The maximal value for \f$t_1\f$
  inline void GetT1extrema(double& t1min_, double& t1max_) const { t1min_=this->fT1min; t1max_=this->fT1max; };
  /// Return the value of the second inner photon's virtuality
  /// \return \f$t_2\f$, the second photon virtuality
  inline double GetT2() const { return this->fT2; };
  /// Return the two limit values for the second photon's virtuality
  /// \param[out] t2min_ The minimal value for \f$t_2\f$
  /// \param[out] t2max_ The maximal value for \f$t_2\f$
  inline void GetT2extrema(double& t2min_, double& t2max_) const { t2min_=this->fT2min; t2max_=this->fT2max; };
  inline double GetS1() const { return this->fS1; };
  inline double GetS2() const { return this->fS2; };
  inline double GetD3() const { return this->_d3; };
  /// Set all the kinematic variables for the outgoing proton remnants, and prepare the hadronisation
  /// \param[in] part_ Particle to "prepare" for the hadronisation to be performed
  void PrepareHadronisation(Particle *part_);
 private:
  /**
   * Calculate energies and momenta of the 1st, 2nd (resp. the "proton-like"
   * and the "electron-like" incoming particles), 3rd (the "proton-like"
   * outgoing particle), 4th (the two-photons central system) and 5th (the
   * "electron-like" outgoing particle) particles in the overall centre-of-mass
   * frame.
   * \brief Energies/momenta computation for the various particles, in the CM
   * system
   * \return Success state of the operation
   */
  bool Orient();
  /**
   * Compute the expression of the matrix element squared for the
   * \f$\gamma\gamma\rightarrow\ell^{+}\ell^{-}\f$ process. It returns the
   * value of the convolution of the form factor or structure functions with
   * the central two-photons matrix element squared.
   * \brief Computes the matrix element squared for the requested process
   * \return Full matrix element for the two-photon production of a pair of
   * spin\f$-\frac{1}{2}-\f$point particles. It is noted as
   * \f[
   *  M = \frac{1}{4bt_1 t_2}\sum_{i=1}^2\sum_{j=1}^2 u_i v_j t_{ij} = \frac{1}{4}\frac{u_1 v_1 t_{11}+u_2 v_1 t_{21}+u_1 v_2 t_{12}+u_2 v_2 t_{22}}{t_1 t_2 b}
   * \f]
   * where \f$b\f$ = \a _bb is defined in \a ComputeWeight as :
   * \f[
   *  b = t_1 t_2+\left(w_{\gamma\gamma}\sin^2{\theta^\text{CM}_6}+4m_\ell\cos^2{\theta^\text{CM}_6}\right) p_g^2
   * \f]
   */
  double PeriPP(int,int);
  /**
   * Describe the kinematics of the process \f$p_1+p_2\to p_3+p_4+p_5\f$ in
   * terms of Lorentz-invariant variables. These variables (along with others)
   * will then be fed into the \a PeriPP method (thus are essential for the
   * evaluation of the full matrix element).
   * \return Success state of the operation
   */
  bool Pickin();
  /// Internal switch for the optimised code version (LPAIR legacy ; unimplemented here)
  int _nOpt;
  
  // --- PICKZZ common block
  
  /// \f$\left|\mathbf p_1\right|\f$, 3-momentum norm of the first proton-like incoming particle
  double _pp1;
  /// \f$E_1\f$, energy of the first proton-like incoming particle
  double _ep1;
  /// \f$m_1\f$, mass of the first proton-like incoming particle
  double _mp1;
  /// \f$m_1^2\f$, squared mass of the first proton-like incoming particle
  double _w1;
  /// \f$\left|\mathbf p_2\right|\f$, 3-momentum norm of the second proton-like incoming particle
  double _pp2;
  /// \f$E_2\f$, energy of the second proton-like incoming particle
  double _ep2;
  /// \f$m_2\f$, mass of the second proton-like incoming particle
  double _mp2;
  /// \f$m_2^2\f$, squared mass of the second proton-like incoming particle
  double _w2;
  /// \f$m_3^2\f$, squared mass of the first proton-like outgoing particle
  double fMX2;
  /// \f$\left|\mathbf p_4\right|\f$, 3-momentum norm of the two-photon central system
  double _pc4;
  /// \f$E_4\f$, energy of the two-photon central system
  double _ec4;
  /// \f$m_4\f$, mass of the two-photon central system
  double _mc4;
  /// \f$m_4^2\f$, squared mass of the two-photon central system
  double _w4;
  /// \f$m_5\f$, mass of the second proton-like outgoing particle
  double fMY2;
  /// \f$m_6^2\f$, squared mass of the first outgoing lepton
  double fMl12;
  /// \f$p_{T,6}\f$, transverse momentum of the first outgoing lepton
  double _pt_l6;
  /// \f$E_6^\mathrm{lab}\f$, energy of the first outgoing lepton, computed in the lab frame
  double _e6lab;
  //// \f$m_7^2\f$, squared mass of the second outgoing lepton
  double fMl22;
  /// \f$p_{T,7}\f$, transverse momentum of the second outgoing lepton
  double _pt_l7;
  /// \f$E_7^\mathrm{lab}\f$, energy of the second outgoing lepton, computed in the lab frame
  double _e7lab;
  /// Energy of the first central photon of momentum \f$t_1\f$
  double _eg1;
  /// 3-momentum of the second central photon of momentum \f$t_1\f$
  double _p3_g1[3];
  /// Energy of the second central photon of momentum \f$t_2\f$
  double _eg2;
  /// 3-momentum of the second central photon of momentum \f$t_2\f$
  double _p3_g2[3];
  /// Total energy provided by the two incoming proton-like particles
  double _etot;
  /// Total momentum provided by the two incoming proton-like particles (along the \f$z\f$-axis)
  double _ptot;
  /// Minimal \f$Q^2\f$ exchange
  double _q2min;
  /// Maximal \f$Q^2\f$ exchange
  double _q2max;
  double _qp2min, _qp2max;
  double _d3;
  
  // --- ACCURA common block
  
  double _acc3;
  double _acc4;
  
  // --- CIVITA common block
  
  double _epsi;
  double _g5, _g6;
  double _a5, _a6;
  double _bb;
  
  // --- DOTP common block
  
  /// \f$p_{12} = \frac{1}{2}\left(s-m_{p_1}^2-m_{p_2}^2\right)\f$
  double _p12;
  /// \f$p_{13} = -\frac{1}{2}\left(t_1-m_{p_1}^2-m_{p_3}^2\right)\f$
  double _p13;
  double _p14, _p15, _p23, _p24, _p25, _p34, _p35, _p45;
  double _p1k2, _p2k1;
  
  // --- DOTPS common block
  
  double _d1dq, _d1dq2, _q1dq, _q1dq2;
  
  // --- EXTRA common block
  
  double fS1, fS2;
  /// Virtuality of the first incoming photon
  double fT1;
  /// Minimal virtuality of the first incoming photon allowed by the kinematics constraints
  double fT1min;
  /// Maximal virtuality of the first incoming photon allowed by the kinematics constraints
  double fT1max;
  /// Virtuality of the second incoming photon
  double fT2;
  /// Minimal virtuality of the second incoming photon allowed by the kinematics constraints
  double fT2min;
  /// Maximal virtuality of the second incoming photon allowed by the kinematics constraints
  double fT2max;
  
  // --- LEVI common block
  
  double _gram;
  double _dd1, _dd2, _dd3, _dd5;
  /**
   * Invariant used to tame divergences in the matrix element computation. It is defined as
   * \f[\Delta = \left(p_1\cdot p_2\right)\left(q_1\cdot q_2\right)-\left(p_1\cdot q_2\right)\left(p_2\cdot q_1\right)\f]
   * with \f$p_i, q_i\f$ the 4-momenta associated to the incoming proton-like particle and to the photon emitted from it.
   */
  double _delta;
  double _g4;
  double _sa1, _sa2;
  
  // --- LEVI common block
  
  /// \f$\delta_1=m_3^2-m_1^2\f$ as defined in Vermaseren's paper
  /// \cite Vermaseren1983347 for the full definition of this quantity
  double _w31;
  double _dw31;
  /// \f$\delta_4=m_5^2-m_2^2\f$ as defined in Vermaseren's paper
  /// \cite Vermaseren1983347 for the full definition of this quantity
  double _w52;
  double _dw52;
  /// \f$\delta_5=m_4^2-t_1\f$ as defined in Vermaseren's paper
  /// \cite Vermaseren1983347 for the full definition of this quantity
  double _dd4;
  /// \f$\delta_2=m_1^2-m_2^2\f$ as defined in Vermaseren's paper
  /// \cite Vermaseren1983347 for the full definition of this quantity
  double _w12;
  
  // --- PICKZZ common block
  
  double _sl1;
  
  // --- VARIAB common block
  
  double _p;
  /// \f$\cos\theta_4\f$ of the two-photons centre-of-mass system
  double _ct4;
  /// \f$\sin\theta_4\f$ of the two-photons centre-of-mass system
  double _st4;
  
  // --- VARIAC common block
  
  double _al4;
  double _be4;
  double _de3, _de5;
  double _pt4;
  
  // --- VARIAD common block
  
  /// Kinematics of the first incoming proton
  Particle::Momentum fP1lab;
  /// Kinematics of the second incoming proton
  Particle::Momentum fP2lab;
  /// Kinematics of the first outgoing proton
  Particle::Momentum fP3lab;
  /// Kinematics of the two-photon system (in the two-proton CM)
  Particle::Momentum fP4lab;
  /// Kinematics of the second outgoing proton
  Particle::Momentum fP5lab;
  /// Kinematics of the first outgoing lepton (in the two-proton CM)
  Particle::Momentum fP6cm;
  /// Kinematics of the second outgoing lepton (in the two-proton CM)
  Particle::Momentum fP7cm;
  double fJacobian;

  double _cotth1, _cotth2;
};

#endif

