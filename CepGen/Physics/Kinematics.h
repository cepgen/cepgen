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
      class Limits : private std::pair<double,double>
      {
        public:
          Limits( double min=invalid_, double max=invalid_ ) : std::pair<double,double>( min, max ) {}
          double lower() const { return first; }
          double upper() const { return second; }
          void in( double low, double up ) { first = low; second = up; }
          double range() const { return second-first; }
          double& lower() { return first; }
          double& upper() { return second; }
          bool hasLower() const { return first != invalid_; }
          bool hasUpper() const { return second != invalid_; }

          friend std::ostream& operator<<( std::ostream&, const Limits& );
        private:
          static constexpr double invalid_ = -999.999;
      };
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
      friend std::ostream& operator<<( std::ostream&, const Cuts& );
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
      friend std::ostream& operator<<( std::ostream&, const ProcessMode& );
  
      /// Dump all the parameters used in this process cross-section computation
      /// or events generation
      void dump( std::ostream& os=std::cout );

      inline void setSqrtS( double sqrts ) { in1p = in2p = sqrts/2; }
      /// First incoming particle's momentum (in \f$\text{GeV}/c\f$)
      double in1p;
      /// Second incoming particle's momentum (in \f$\text{GeV}/c\f$)
      double in2p;
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
      ProcessMode mode;
      StructureFunctions remnant_mode;
      /// Sets of cuts to apply on the final phase space
      Cuts cuts_mode;
      /// Limits on the transverse momentum of the single outgoing particles
      Limits pt_single_central;
      /// Limits on the energy of the central two-photons system
      Limits e_single_central;
      /// Limits on the pseudo-rapidity (\f$\eta\f$) of the outgoing particles
      Limits eta_single_central;
      /// Limits on the mass of the central system
      Limits mass_central;
      /// Limits on the mass (in GeV/c\f${}^\mathrm{2}\f$) of the outgoing proton remnant(s)
      Limits mass_remnants;
      /// Limits on the value of \f$Q^2\f$
      Limits q2;
      /// Limits on \f$s\f$ on which the cross section is integrated
      Limits w;
      /// Limits on the difference in outgoing particles' transverse momentum
      Limits pt_diff_central;
      /// Limits on the transverse component of the energy transfer
      Limits qt;
  };
}

#endif

