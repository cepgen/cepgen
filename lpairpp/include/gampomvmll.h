#ifndef _GAMPOMVMLL_H
#define _GAMPOMVMLL_H

#include <algorithm>

#include "process.h"

/**
 * @brief Computes the matrix element for a CE \f$\gamma\mathbb{P}\rightarrow J/\psi,\Upsilon\rightarrow\ell^+\ell^-\f$ process
 */
class GamPomVMLL : public Process
{
 public:
  GamPomVMLL();
  ~GamPomVMLL();
  bool SetIncomingParticles(Particle, Particle);
  bool SetOutgoingParticles(int, int);
  void FillKinematics(bool);
  void SetKinematics(Kinematics);
  void ComputeCMenergy();
  double ComputeMX(double x_, double outmass_, double* dw_);
  double ComputeWeight();
  void StoreEvent(std::ofstream*,double);
  void PrepareHadronisation(Particle *part_);
 private:
  void GenGam();
  /**
   * Generate one event with unweighted photon & electron
   * * according to WWA :
   *   - transversal photonspectrum. \f$Q^2\rightarrow 0\f$:
   *
   *     \f$P(y,Q^2)=\frac{\alpha}{2\pi}\frac{1}{Q^2y}\left(2(1-y)\left(1-\frac{Q^2_\mathrm{min}}{Q^2}\right)+y^2\right)\f$
   *
   *   - longitudinal photonspectrum. \f$Q^2\rightarrow 0\f$:
   *
   *     \f$P(y,Q^2)=\frac{\alpha}{2\pi}\frac{1}{Q^2y}\left(2(1-y)\right)\f$
   *
   * * full transversal photonspectrum given by:
   *   - ABT, I. & J.R. SMITH (1992): MC upgrades to study untagged events. - H1-10/92-249.
   *   - SMITH, J.R. (1992): An experimentalist's guide to photon flux calculations. - H1-12/92-259
   *   - SMITH, J.R. & B.D. BUROW (1994): Photon fluxes with beam mass effects and polarizations. - H1-01/94-338.
   * * full transversal and longitudinal spectrum by ABT&SMITH
   *   - calculate integrated factor over the spectrum:
   *     kinematical bounds : \f$\left[Y_{\mathrm{min}},  Y_{\mathrm{max}}\right] (W_{\mathrm{min}})\f$, 
   *                       \f$\left[Q^2_{\mathrm{min}}, Q^2_{\mathrm{max}}\right] (Q^2_{\mathrm{cutoff}})\f$
   * @param igammd_ Photon generation mode:
   *    - 1: WWA/EPA approximation (including electron-mass effect and longitudinal flux). **Recommended**
   *    - 2: Transverse spectrum
   *    - 3: Transverse & longitudinal spectrum
   *    - 4: as 3, but flux in proton rest frame
   */
  void GEPhot(int igammd_=1);
  double PXMass(double mmin_, double mmax_);
  /**
   * Generate hadronic mass between @a mmin_ and @a mmax_ for VM vertex
   * @param mmin_ Minimal allowed mass
   * @param mmax_ Maximal allowed mass
   * @return Hadronic mass in GeV
   * @author Benno List
   * @date 14 jan 1992
   */
  double VXMass(double mmin_, double mmax_);
  void FragGl();

  /**
   * Let the generated vector meson decay
   * @author Benno List
   * @date 25 jan 1993
   */
  void DecVM();

  std::string _name;
  int ifragp, ifragv;
  /**
   * PDG code for produced vector meson (should have \f$J^{PC}=1^{--}\f$)
   * Possible values:
   *  - 113 : \f$\rho\f$
   *  - 223 : \f$\omega\f$
   *  - 333 : \f$\phi\f$
   *  - 443 : \f$J/\psi\f$
   *  - 20443 : \f$\psi'\f$
   *  - 553 : \f$\Upsilon(1s)\f$
   *  - 20553 : \f$\Upsilon(2s)\f$
   *  - 30553 : \f$\Upsilon(3s)\f$
   *  - 40113 : \f$\rho(1450)\rightarrow\pi^+\pi^-\rho^0\f$
   *  - 10333 : \f$\phi(1680)\rightarrow K\bar K\f$
   *  - *22 : diffr. gamma dissoc. (special value)*
   * @brief PDG code for produced vector meson
   */
  int itypvm;

  double pe, dme, pp, dmp; //FIXME

  //// Parameters for the pomeron
  /**
   * @brief Intercept of pomeron trajectory minus 1
   * @note Controls rise of \f$\sigma_{\gamma p}\f$ with \f$W\f$
   */
  double _epsilw;
  /**
   * @brief Intercept of pomeron trajectory minus 1
   * @note Controls \f$M_X\f$ spectrum
   */
  double _epsilm;
  /**
   * @brief Slope \f$\alpha'\f$ of pomeron trajectory in \f$\mathrm{GeV}^{-2}\f$
   * @note Controls shrinkage of \f$b\f$ slope
   */
  double _alpha1;


  double _s, _ecm;

  //// Photon generator mode
  /**
   * @brief Minimal CM energy of \f$\gamma p\f$ system
   */
  double _wmin;
  /**
   * @brief Maximal CM energy of \f$\gamma p\f$ system
   * @note If too low, it is set to \f$\sqrt s\f$
   */
  double _wmax;
  /**
   * @brief Minimal \f$Q^2\f$ of photon in \f$\mathrm{GeV}^2\f$
   */
  double _q2min;
  /**
   * @brief Maximal \f$Q^2\f$ of photon in \f$\mathrm{GeV}^2\f$
   */
  double _q2max;
  /**
   * @brief Minimal value of scaling variable \f$y\f$
   */
  double _ymin;
  /**
   * @brief Maximal value of scaling variable \f$y\f$
   */
  double _ymax;

  //// Parameters for t spectrum
  /**
   * Slope parameter \f$b\f$ of \f$t\f$ distribution in \f$\mathrm{GeV}^{-2}\f$
   * at CM energy @a _wb0 and (for diffractive dissociation) mass @a _amxb0
   * @note Must be positive!
   */
  double _b0;
  /**
   * @brief CM energy of \f$\gamma p\f$ system at which @a _b0 was measured,
   * in GeV
   */
  double _wb0;
  /**
   * @brief Mass of diffractively dissociating hadronic system for which
   * @a _b0 was measured
   * @note For @a _amxb0=0.0, @a _amxb0 is set according to production mode.
   *  Value is not meaningful for elastic VM production
   */
  double _amxb0;
  /**
   * Power law exponent.
   *  - For @a _anexp\f$=0\f$ (default), a pure exponential spectrum is
   *    generated according to  (taking t<0)
   *  \f$\frac{\mathrm d\sigma}{\mathrm dt} = e^{bt}\f$
   *  - For @a _anexp\f$>1\f$, an interpolated spectrum is generated
   *    according to
   *  \f$\frac{\mathrm d\sigma}{\mathrm dt} = \exp\left[-n\ln{\left(-\frac{bt}{n}+1\right)}\right] = \left(-\frac{bt}{n}+1\right)^{-n}\f$
   *    with \f$n=\f$ @a _anexp
   *    - Limit for small \f$bt\f$: \f$\exp\left(bt+ct^2\right)\f$ with \f$c=b^2/2n\f$
   *    - Limit for large \f$bt\gg n\f$: \f$t^{-n}\f$
   */
  double _anexp;

  /**
   * @brief \f$\gamma p\f$ CM energy at which SIGGP was measured
   */
  double _wsig0;
  double _w2;
  /**
   * @brief Absolute of square-momentum of virtual photon
   */
  double _q2;

  /**
   * @brief CM momentum of outgoing particles
   */
  double _pcm3;
  double _pcmvm[3];

  bool _genmxt_begin;

  bool _gengam_first;

  bool _gephot_first;
  double _gephot_pel[5];
  double _gephot_ppr[5];
  double _gephot_pph[5];
  double _gephot_ppe[5];
  double _gephot_q2;
  int _gephot_heli;

  bool _fraggl_begin;

  double _ppcms8[1000][5]; // in cprod : PPCMS8:  5-vectors of particles in the gamma-p CMS in double prec.

  /**
   * @brief Mass of generated vector meson
   */
  double _dmvm;
  /**
   * @brief Width of generated vector meson
   */
  double _dwvm;
};

#endif
