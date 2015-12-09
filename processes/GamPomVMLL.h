#ifndef GamPomVMLL_h
#define GamPomVMLL_h

#include <algorithm>

#include "../include/utils.h"
#include "../include/GenericProcess.h"

#define IBE 1
#define ISCE 3
#define IGAM 41
#define IVVM 411

#define IBP 2
#define IDIFP 5
#define IPOM 42

#define IDIFV 43
#define IGLUE 431
#define IVM 4

#define OL1 6
#define OL2 7

/**
 * @brief Computes the matrix element for a CE \f$\gamma\mathbb{P}\rightarrow \rho,\omega,\phi,J/\psi,\Upsilon,\ldots\rightarrow\ell^+\ell^-\f$ process
 */
class GamPomVMLL : public GenericProcess
{
 public:
  GamPomVMLL();
  ~GamPomVMLL();
  void AddEventContent();
  //bool SetIncomingParticles(Particle, Particle);
  //bool SetOutgoingParticles(int, int);
  void FillKinematics(bool);
  //void SetKinematics(Kinematics);
  //void ComputeCMenergy();
  double ComputeWeight();
  //void StoreEvent(std::ofstream*,double);
  //void PrepareHadronisation(Particle *part_);
 private:
 /**
  * Set up the generator for event generation
  * @author Benno List
  * @date 31 Jan 1993 (INIGEN)
  * @date 14 Jun 1995 (GDIBEG)
  * @note Some of the parameters are set/calculated/changed during this step
  */
  void GDIBeg();
  void GDIEvt();
  /**
   * Generate a diffractive vector meson production event
   */
  void GenEvtDi();
  /**
   * Take 5-vectors of colliding \f$e\f$ and \f$p\f$ and generate a virtual photon, the momentum transfer of the pomeron, and the diffractive masses at the \f$p\f$ and VM vertices
   */
  void GenGam();
  double OneEvent();

  inline void GenBEl() {};
  inline void GenBPr() {};
  /**
   * Generate \f$m_X^p\f$, \f$m_X^\text{VM}\f$ and \f$t\f$ and determine if the generated combination is kinematically allowed
   * @author Benno List
   * @date 18 Jul 1993 (MXTGEN)
   * @date 10 May 1994 (GENMXT)
   */
  double GenMXT(double* wght);
  /**
   * Take 5-vectors of colliding \f$\gamma\f$ and \f$p\f$ and generate a diffractive state
   * @param[in] tpom_ Momentum transfer of the pomeron
   * @param[in] yhat_ \f$\hat y=\sin^2(\theta^\ast/2)\f$ where \f$\theta^\ast\f$ is the scattering angle in the \f$\gamma p\f$ CMS
   */
  void GenDif();
  /**
   * Generate one event with unweighted photon & electron
   * * according to WWA :
   *   - transversal photonspectrum. \f$Q^2\rightarrow 0\f$:
   *     \f[P(y,Q^2)=\frac{\alpha}{2\pi}\frac{1}{Q^2y}\left(2(1-y)\left(1-\frac{Q^2_\text{min}}{Q^2}\right)+y^2\right)\f]
   *   - longitudinal photonspectrum. \f$Q^2\rightarrow 0\f$:
   *     \f[P(y,Q^2)=\frac{\alpha}{2\pi}\frac{1}{Q^2y}\left(2(1-y)\right)\f]
   *
   * * full transversal photonspectrum given by @cite AbtSmith:1992, @cite Smith:1992, @cite SmithBurow:1994
   * * full transversal and longitudinal spectrum by @cite AbtSmith:1992
   *   - calculate integrated factor over the spectrum:
   *     kinematical bounds : \f$\left[Y_{\text{min}},  Y_{\text{max}}\right] (W_{\text{min}})\f$, 
   *                       \f$\left[Q^2_{\text{min}}, Q^2_{\text{max}}\right] (Q^2_{\text{cutoff}})\f$
   * @param[out] q2_ Virtuality of photon (positive!): \f$Q^2 = -q^2\f$
   * @param[out] heli_ Photon helicity: 0: longitudinal, 1, -1: transverse polarization
   * @author T. Jansen
   * @date 6 Apr 1993
   */
  void GEPhot(double* q2_, int* heli_);
  double PXMass(double mmin_, double mmax_);
  /**
   * Generate hadronic mass between @a mmin_ and @a mmax_ for VM vertex
   * @param[in] mmin_ Minimal allowed mass
   * @param[in] mmax_ Maximal allowed mass
   * @return Hadronic mass in GeV
   * @author Benno List
   * @date 14 Jan 1992
   */
  double VXMass(double mmin_, double mmax_);
  void FragGl();

  /**
   * Let the generated vector meson decay
   * @author Benno List
   * @date 25 Jan 1993
   */
  void DecVM();
  /**
   * Generate photon with energy between @a emin_ and electron energy and \f$Q^2\f$ less than @a q2max_, calculate 5-vector of scattered electron.
   * @param[in] emin_ Minimal allowed energy
   * @param[in] q2max_ Maximal allowed \f$Q^2\f$
   * @param[in] pel_ Incoming electron kinematics
   * @param[out] phot_ Photon kinematics (mass is set to \f$-\sqrt Q^2\f$)
   * @param[out] ele_ Outgoing electron kinematics
   * @param[out] q2_ Absolute Photon momentum squared (positive!)
   * @author Benno List
   * @date 22 Jan 1993
   * @note Up to now only real photons with 1/k spectrum
   * @note Sets @a _ihel, @a _ftrans, @a _epsil
   */
  void GenPhot(Particle* phot_, Particle* ele_, double *q2_, Particle pel_, double emin_, double q2max_);
  /**
   * Generate photon with fixed energy @a _egamma, and calculate scattered electron kinematics.
   * @param[in] pel_ Incoming electron kinematics
   * @param[in] egamma_ Photon energy
   * @param[out] phot_ Photon kinematics (mass is set to \f$-\sqrt Q^2\f$)
   * @param[out] ele_ Outgoing electron kinematics
   * @param[out] q2_ Absolute Photon momentum squared (positive!)
   * @author Benno List
   * @date 28 Apr 1993
   */
  void FixPhot(Particle* phot_, Particle* ele_, double *q2_, Particle pel_, double egamma_);
  /**
   * Calculate relative photon luminosity for photon flux produced by @a GEPhot, weighted by VM propagator and cross section
   * @return _genmxt_f Total VM flux, relative to \f$e\f$ flux
   * @return _genmxt_df Error of @a _genmxt_f
   * @return _genmxt_ft Relative VM flux for transverse VMs
   * @return _genmxt_dft Error of @a _genmxt_ft
   * @return _genmxt_fl Relative VM flux for longitudinal VMs
   * @return _genmxt_dfl Error of @a _genmxt_fl
   * @author T. Jansen 
   * @date 07 Apr 1993
   */
  void VMFlux();

  int _event_heli;
  double _event_egammin;
  double _event_smax;
  double _event_propmx;

  /**
   * Minimal \f$\cos\theta\f$ of scattered electron
   */
  double _cthelb;
  /**
   * Minimal energy of scattered electron in GeV
   */
  double _eelmin;

  /**
   * Fragmentation mode for diffractive proton state.
   * Possible values:
   *  - 0 : Elastic scattering of proton
   *  - 1 : Fragmentation by JETSET 7.3 with gluon emission (called DIFFVMg, contributed by Leszek Adamczyk)
   *  - 1 : Fragmentation by JETSET 7.3
   *  - 2 : Isotropic phase space decay into nucleon and pions
   *  - 12212 : Elastic \f$N(1440)^+\f$ production at \f$p\f$ vertex (for other \f$N^\ast\f$ states, insert respective PDG code)
   */
  int ifragp;
  /**
   * Minimal energy released in decay of diffractive proton state, in GeV
   * If value is too small (smaller than \f$m_{\pi^0}\f$), sets value to \f$0.236 = m_n+m_{\pi^0}+0.10-m_p\f$
   * @note Value is only meaningful for @a ifragp \f$= 1$ or $2\f$.
   */
  double deminp;
  /**
   * Fragmentation mode for diffractive vector meson state.
   * Possible values:
   *  - 0 : Elastic vector meson production
   *  - 1 : Fragmentation by JETSET 7.3
   *  - 2 : Isotropic phase space decay into VM+pions
   *  - 995 : diffractive pomeron-VM scattering (glueball production)
   *    see P.E.SCHLEIN (1994): Phys. Lett. B332, 136-140.
   */
  Particle::ParticleCode ifragv;
  /**
   * Minimal mass of diffractive VM state.
   * If value is too small (smaller than \f$m_{\pi^0}\f$), sets value to \f$m_{\text{VM}}+\f$ some offset
   * @note Value is only meaningful for @a ifragv \f$= 1\f$ or \f$2\f$ or \f$955\f$.
   */
  double amassv;

  ////
  /**
   * Type of vector meson (should have \f$J^{PC}=1^{--}\f$) to produce, and decay mode
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
   * @brief Type of vector meson to produce and its decay channel
   */
  Particle::ParticleCode itypvm;
  /**
   * @brief Index of diffractive \f$q\bar q\f$ states
   */
  int idifv;
  /**
   * @brief Index of virtual vector meson
   */
  int ivvm;
  /**
   * @brief Index of pomeron photon
   */
  int ipom;
  /**
   * @brief Index of vector meson
   */
  int ivm;

  //// VMD model parameters
  /**
   * Parameter for \f$Q^2\f$-dependence of cross section in GeV:
   * \f[\sigma(Q^2) = \frac{\sigma(0)}{\left(1+\frac{Q^2}{\Lambda^2}\right)^{\epsilon_\text{prop}}}\f]
   */
  double _lambda;
  /**
   * Propagator term exponent \f$\epsilon_\text{prop}\f$ (see @a _lambda)
   */
  double _eprop;
  /**
   * Parameter for \f$Q^2\f$-dependence of \f$\sigma_L/\sigma_T\f$:
   * \f[
   *  \frac{\sigma_L(Q^2)}{\sigma_T(Q^2)} = \frac{\xi Q^2/m^2}{1+\xi\chi Q^2/m^2}
   * \f]
   * with
   * - \f$\frac{\sigma_L}{\sigma_T}\rightarrow\xi\frac{Q^2}{m^2}\f$ for low \f$Q^2\f$
   * - \f$\frac{\sigma_L}{\sigma_T}\rightarrow\frac{1}{\chi}\f$ for high \f$Q^2\f$
   *
   * @a _xi is assumed to be less than 4 (more precisely, it is assumed that \f$\sigma_L(Q^2)\f$ is always less than \f$\sigma_T(0)\f$).
   */
  double _xi;
  /**
   * See @a _xi. \f$\chi\f$ is a purely phenomenological parameter with no theoretical justification!
   */
  double _chi;

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
   * @brief Slope \f$\alpha'\f$ of pomeron trajectory in \f$\text{GeV}^\text{-2}\f$
   * @note Controls shrinkage of \f$b\f$ slope
   */
  double _alpha1;
  double _alph1m;

  //// Photon generator mode
  /**
   * Photon generator mode.
   * Possible values:
   *    - -1: Fixed photon energy @a _egamma
   *    -  0: \f$\frac{1}{k}\f$ spectrum
   *    - 1: WWA/EPA approximation (including electron-mass effect and longitudinal flux). **Recommended**
   *    - 2: Transverse spectrum _a la_ @cite AbtSmith:1992
   *    - 3: Transverse & longitudinal spectrum _a la_ @cite AbtSmith:1992
   *    - 4: as 3, but flux in proton rest frame
   */
  int _igammd;
  /**
   * Energy of photon in GeV for @a _igammd \f$= -1\f$
   */
  double _egamma;
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
   * @brief Minimal \f$Q^2\f$ of photon in \f$\text{GeV}^\text{2}\f$
   */
  double _q2min;
  /**
   * @brief Maximal \f$Q^2\f$ of photon in \f$\text{GeV}^\text{2}\f$
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
   * Slope parameter \f$b\f$ of \f$t\f$ distribution in \f$\text{GeV}^\text{-2}\f$
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
   *    generated according to  (taking \f$t<0\f$)
   *  \f[\frac{\text d\sigma}{\text dt} = e^{bt}\f]
   *  - For @a _anexp\f$>1\f$, an interpolated spectrum is generated
   *    according to
   *  \f[\frac{\text d\sigma}{\text dt} = \exp\left[-n\ln{\left(-\frac{bt}{n}+1\right)}\right] = \left(-\frac{bt}{n}+1\right)^{-n}\f]
   *    with \f$n=\f$ @a _anexp
   *    - Limit for small \f$bt\f$: \f[\exp\left(bt+ct^2\right)\f] with \f$c=b^2/2n\f$
   *    - Limit for large \f$bt\gg n\f$: \f$t^{-n}\f$
   */
  double _anexp;

  /**
   * @brief \f$\gamma p\f$ CM energy at which SIGGP was measured
   */
  double _wsig0;
  /**
   * Branching ratio of the chosen decay channel.
   * Useful values:
   *  - 1      for @a itypel \f$= 0\f$
   *  - 0.99   for \f$\rho^0\rightarrow \pi^+ \pi^-\f$
   *  - 0.0221 for \f$\omega\rightarrow \pi^+ \pi^-\f$
   *  - 0.491  for \f$\phi\rightarrow K^+ K^-\f$
   *  - 0.344  for \f$\phi\rightarrow K^0_L K^0_S\f$
   *  - 0.0598 for \f$J/\psi\rightarrow e^+ e^-, \mu^+ \mu^-\f$
   *  - 0.0425 for \f$\psi'\rightarrow \ell^+ \ell^- X\f$ (including cascade decays)
   *  - 0.025  for \f$\Upsilon(1s)\rightarrow \ell^+ \ell^-\f$
   *  - 0.02   for \f$\Upsilon(2s)\rightarrow \ell^+ \ell^- X\f$ (including cascade decays)
   *  - 0.0217 for \f$\Upsilon(3s)\rightarrow \ell^+ \ell^- X\f$ (including cascade decays)
   */
  double _br;
  double _gengam_w2;
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
  double _genmxt_b;

  bool _gengam_first;
  double _gengam_yhat;
  double _gengam_t;

  //static bool _gephot_first; // FIXME FIXME no static variables should be allowed here !
  bool _gephot_first;
  double _gephot_pel[5];
  double _gephot_ppr[5];
  double _gephot_pph[5];
  double _gephot_ppe[5];
  int _gephot_heli;

  bool _fraggl_begin;

  //// Common block /PHOTINT/
  double _photint_swei, _photint_swei2, _photint_sweit, _photint_sweit2, _photint_sweil, _photint_sweil2;

  double _ppcms8[1000][5]; // in cprod : PPCMS8:  5-vectors of particles in the gamma-p CMS in double prec.

  /**
   * @brief Mass of generated vector meson
   */
  double _dmvm;
  /**
   * @brief Width of generated vector meson
   */
  double _dwvm;

  /**
   * @brief Mass at the proton-pomeron vertex
   */
  double _genmxt_dmxp;
  /**
   * @brief Mass at the vector meson-pomeron vertex
   */
  double _genmxt_dmxv;
  double _genmxt_bmin;

  double _dme;
  double _dmp;
  double _dmpi;
  double _dmpi0;
  double _dmn;
  double _dml;
  double _dmnst;
  double _dwnst;
  double _pz1;
  double _e1;
  double _pz2;
  double _e2;
  
  double _vmflux_f;
  double _vmflux_df;
  double _vmflux_fl;
  double _vmflux_dfl;
  double _vmflux_ft;
  double _vmflux_dft;
  int _iacct;
  int _iaccl;
  int _isum;
  int _igen;
  int _igent;
  int _igenl;
  double _qsumt;
  double _qsuml;
  double _dsumt;
  double _dsuml;

};

#endif
