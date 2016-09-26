#ifndef Kinematics_h
#define Kinematics_h

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>

#include "utils.h"

/// List of kinematic cuts to apply on the central and outgoing phase space.
class Kinematics
{
 public:
  Kinematics();
  ~Kinematics();

  /**
   * @brief Set of cuts to apply on the central system
   * - 0 - No cuts at all (for the total cross section)
   * - 1 - Vermaserens' hypothetical detector cuts : for both leptons,
   *   + \f$\frac{|p_z|}{|\mathbf p|}\leq\f$ 0.75 and \f$p_T\geq 1~\text{GeV}/c\f$,
   *   or
   *   + 0.75 \f$<\frac{|p_z|}{|\mathbf p|}\leq\f$ 0.95 and \f$p_z> 1~\text{GeV}/c\f$,
   * - 2 - Cuts on both the outgoing leptons, according to the provided cuts parameters
   * - 3 - Cuts on at least one outgoing lepton, according to the provided cut parameters
   */
  enum Cuts { NoCuts = 0, VermaserenCuts = 1, BothLeptons = 2, OneLepton = 3 };
  /// Type of outgoing process kinematics to be considered (elastic/dissociative final states)
  enum ProcessMode {
    ElectronProton = 0,
    ElasticElastic = 1,
    ElasticInelastic = 2,
    InelasticElastic = 3,
    InelasticInelastic = 4
  };
  /// Human-readable format of a process mode (elastic/dissociative parts)
  friend std::ostream& operator<<(std::ostream& os, const Kinematics::ProcessMode& pm);
  
  /// Dump all the parameters used in this process cross-section computation
  /// or events generation
  void Dump();
  /**
   * Type of kinematics to consider for the process. Can either be :
   *  * 0 for the electron-proton elastic case
   *  * 1 for the proton-proton elastic case
   *  * 2 for the proton-proton single-dissociative (or inelastic-elastic) case
   *  * 3 for the proton-proton single-dissociative (or elastic-inelastic) case
   *  * 4 for the proton-proton double-dissociative case
   * @brief Type of kinematics to consider for the phase space
   */
  ProcessMode kinematics;
  /// Sets of cuts to apply on the final phase space
  Cuts mode;
  /// Minimal transverse momentum of the single outgoing leptons
  double ptmin;
  /// Maximal transverse momentum of the single outgoing leptons
  double ptmax;
  /// Minimal energy of the central two-photons system
  double emin;
  /// Maximal energy of the central two-photons system
  double emax;
  /// Minimal rapidity (\f$\eta\f$) of the outgoing lepton
  double etamin;
  /// Maximal rapidity (\f$\eta\f$) of the outgoing lepton
  double etamax;
  /// Minimal mass (in GeV/c\f${}^\mathrm{2}\f$) of the outgoing proton remnant(s)
  double mxmin;
  /// Maximal mass (in GeV/c\f${}^\mathrm{2}\f$) of the outgoing proton remnant(s)
  double mxmax;
  /// Minimal value of \f$Q^2\f$
  double q2min;
  /// Maximal value of \f$Q^2\f$
  double q2max;
  /// Minimal \f$s\f$ on which the cross section is integrated
  double wmin;
  /// Maximal \f$s\f$ on which the cross section is integrated. If negative,
  /// the maximal energy available to the system (hence, \f$s=(\sqrt{s})^{2}\f$)
  /// is provided.
  double wmax;
  /// Minimal difference in outgoing particles' transverse momentum
  double ptdiffmin;
  /// Maximal difference in outgoing particles' transverse momentum
  double ptdiffmax;
  double qtmin;
  double qtmax;
};

#endif

