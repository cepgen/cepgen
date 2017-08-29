#ifndef CepGen_Physics_Kinematics_h
#define CepGen_Physics_Kinematics_h

#include <iomanip>
#include <algorithm>
#include <string>

#include "CepGen/Core/utils.h"
#include "StructureFunctions.h"
#include "Particle.h"
#include "Cuts.h"

using std::cout;

namespace CepGen
{
  /// List of kinematic cuts to apply on the central and outgoing phase space.
  class Kinematics
  {
    public:
      /// Validity interval for a variable
      class Limits : private std::pair<double,double>
      {
        public:
          Limits( double min=invalid_, double max=invalid_ ) : std::pair<double,double>( min, max ) {}

          /// Lower limit to apply on the variable
          double min() const { return first; }
          /// Lower limit to apply on the variable
          double& min() { return first; }
          /// Upper limit to apply on the variable
          double max() const { return second; }
          /// Upper limit to apply on the variable
          double& max() { return second; }
          /// Specify the lower and upper limits on the variable
          void in( double low, double up ) { first = low; second = up; }
          /// Full variable range allowed
          double range() const { return ( !hasMin() || !hasMax() ) ? 0. : second-first; }
          /// Have a lower limit?
          bool hasMin() const { return first != invalid_; }
          /// Have an upper limit?
          bool hasMax() const { return second != invalid_; }

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
      enum CutsMode { NoCuts = 0, BothParticles = 2, OneParticle = 3 };
      /// Human-readable format of a cuts mode
      friend std::ostream& operator<<( std::ostream&, const CutsMode& );
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
      void dump( std::ostream& os=std::cout ) const;

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
      CutsMode cuts_mode;
      std::map<Cuts::Central, Limits> central_cuts;
      std::map<Cuts::Remnants, Limits> remnant_cuts;
      std::map<Cuts::InitialState, Limits> initial_cuts;
  };
}

#endif

