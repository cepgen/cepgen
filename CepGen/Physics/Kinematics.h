#ifndef CepGen_Physics_Kinematics_h
#define CepGen_Physics_Kinematics_h

#include <iomanip>
#include <algorithm>
#include <string>

#include "CepGen/Core/utils.h"
#include "StructureFunctions.h"
#include "Particle.h"

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

      inline void setSqrtS( double sqrts ) { in1p = in2p = sqrts/2; }
      /// First incoming particle's momentum (in \f$\text{GeV}/c\f$)
      float in1p;
      /// Second incoming particle's momentum (in \f$\text{GeV}/c\f$)
      float in2p;
      /// First beam/primary particle's PDG identifier
      Particle::ParticleCode in1pdg;
      /// Second beam/primary particle's PDG identifier
      Particle::ParticleCode in2pdg;
      /// PDG id of the outgoing central particles
      Particle::ParticleCode pair;

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
      Cuts cuts_mode;
      /// Minimal transverse momentum of the single outgoing particles
      double pt_min;
      /// Maximal transverse momentum of the single outgoing particles
      double pt_max;
      /// Minimal energy of the central two-photons system
      double e_min;
      /// Maximal energy of the central two-photons system
      double e_max;
      /// Minimal rapidity (\f$\eta\f$) of the outgoing particles
      double eta_min;
      /// Maximal rapidity (\f$\eta\f$) of the outgoing particles
      double eta_max;
      /// Minimal mass of the central system
      double mass_min;
      /// Minimal mass of the central system
      double mass_max;
      /// Minimal mass (in GeV/c\f${}^\mathrm{2}\f$) of the outgoing proton remnant(s)
      double mx_min;
      /// Maximal mass (in GeV/c\f${}^\mathrm{2}\f$) of the outgoing proton remnant(s)
      double mx_max;
      /// Minimal value of \f$Q^2\f$
      double q2_min;
      /// Maximal value of \f$Q^2\f$
      double q2_max;
      /// Minimal \f$s\f$ on which the cross section is integrated
      double w_min;
      /// Maximal \f$s\f$ on which the cross section is integrated. If negative,
      /// the maximal energy available to the system (hence, \f$s=(\sqrt{s})^{2}\f$)
      /// is provided.
      double w_max;
      /// Minimal difference in outgoing particles' transverse momentum
      double ptdiff_min;
      /// Maximal difference in outgoing particles' transverse momentum
      double ptdiff_max;
      /// Minimal transverse component of the energy transfer
      double qt_min;
      /// Maximal transverse component of the energy transfer
      double qt_max;
  };
}

#endif

