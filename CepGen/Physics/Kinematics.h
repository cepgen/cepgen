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
  /// List of kinematic constraints to apply on the process phase space.
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

      /// Set of cuts to apply on the central system
      enum CutsMode {
        NoCuts = 0,       ///< No cuts at all (for the total cross section)
        AllParticles = 2, ///< Cuts on all the outgoing central particles
        OneParticle = 3   ///< Cuts on at least one outgoing central particle
      };
      /// Human-readable format of a cuts mode
      friend std::ostream& operator<<( std::ostream&, const CutsMode& );
      /// Type of kinematics to consider for the process
      enum ProcessMode {
        ElectronProton = 0,     ///< electron-proton elastic case
        ElasticElastic = 1,     ///< proton-proton elastic case
        ElasticInelastic = 2,   ///< proton-proton single-dissociative (or inelastic-elastic) case
        InelasticElastic = 3,   ///< proton-proton single-dissociative (or elastic-inelastic) case
        InelasticInelastic = 4, ///< proton-proton double-dissociative case
        ProtonElectron,
        ElectronElectron
      };
      /// Human-readable format of a process mode (elastic/dissociative parts)
      friend std::ostream& operator<<( std::ostream&, const ProcessMode& );
  
      /// Dump all the parameters used in this process cross-section computation
      /// or events generation
      void dump( std::ostream& os=std::cout ) const;

      /// Incoming particles' momentum (in \f$\text{GeV}/c\f$)
      std::pair<double,double> inp;
      inline void setSqrtS( double sqrts ) { inp = { sqrts*0.5, sqrts*0.5 }; }
      /// Beam/primary particle's PDG identifier
      std::pair<Particle::ParticleCode,Particle::ParticleCode> inpdg;
      /// PDG id of the outgoing central particles
      std::vector<Particle::ParticleCode> central_system;

      /// Type of kinematics to consider for the phase space
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

