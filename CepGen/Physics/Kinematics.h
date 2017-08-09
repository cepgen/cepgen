#ifndef Kinematics_h
#define Kinematics_h

#include <iomanip>
#include <algorithm>
#include <string>

#include "CepGen/Core/utils.h"
#include "StructureFunctions.h"

using std::cout;

namespace CepGen
{
  /// List of kinematic cuts to apply on the central and outgoing phase space.
  class Kinematics
  {
    public:
      Kinematics();
      ~Kinematics();

      /**
       * \brief Set of cuts to apply on the central system
       * - 0 - No cuts at all (for the total cross section)
       * - 2 - Cuts on both the outgoing central particles, according to the provided cuts parameters
       * - 3 - Cuts on at least one outgoing central particle, according to the provided cut parameters
       */
      enum Cuts { NoCuts = 0, BothParticles = 2, OneParticle = 3 };
      /// Human-readable format of a cuts mode
      friend std::ostream& operator<<( std::ostream&, const Kinematics::Cuts& );
      /// Type of outgoing process kinematics to be considered (elastic/dissociative final states)
      enum ProcessMode {
        ElectronProton = 0,
        ElasticElastic = 1,
        ElasticInelastic = 2,
        InelasticElastic = 3,
        InelasticInelastic = 4,
        ProtonElectron,
        ElectronElectron
      };
      /// Human-readable format of a process mode (elastic/dissociative parts)
      friend std::ostream& operator<<( std::ostream&, const Kinematics::ProcessMode& );
  
      /// Dump all the parameters used in this process cross-section computation
      /// or events generation
      void dump( std::ostream& os=std::cout );
      /**
       * Type of kinematics to consider for the process. Can either be :
       *  * 0 for the electron-proton elastic case
       *  * 1 for the proton-proton elastic case
       *  * 2 for the proton-proton single-dissociative (or inelastic-elastic) case
       *  * 3 for the proton-proton single-dissociative (or elastic-inelastic) case
       *  * 4 for the proton-proton double-dissociative case
       * \brief Type of kinematics to consider for the phase space
       */
      ProcessMode kinematics;
      StructureFunctions remnant_mode;
      /// Sets of cuts to apply on the final phase space
      Cuts mode;
      /// Minimal transverse momentum of the single outgoing particles
      double ptmin;
      /// Maximal transverse momentum of the single outgoing particles
      double ptmax;
      /// Minimal energy of the central two-photons system
      double emin;
      /// Maximal energy of the central two-photons system
      double emax;
      /// Minimal rapidity (\f$\eta\f$) of the outgoing particles
      double etamin;
      /// Maximal rapidity (\f$\eta\f$) of the outgoing particles
      double etamax;
      /// Minimal mass of the central system
      double massmin;
      /// Minimal mass of the central system
      double massmax;
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
      /// Minimal transverse component of the energy transfer
      double qtmin;
      /// Maximal transverse component of the energy transfer
      double qtmax;
  };
}

#endif

