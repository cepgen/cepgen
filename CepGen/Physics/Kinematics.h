#ifndef CepGen_Physics_Kinematics_h
#define CepGen_Physics_Kinematics_h

#include <iomanip>
#include <algorithm>

#include "CepGen/Core/Logger.h"
#include "CepGen/StructureFunctions/StructureFunctions.h"
#include "CepGen/Physics/ParticleProperties.h"

#include "Cuts.h"

#include <vector>
#include <map>

using std::cout;
using std::string;

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
          /// Define lower and upper limits on a quantity
          Limits( double min = invalid_, double max = invalid_ ) :
            std::pair<double,double>( min, max ) {}

          /// Lower limit to apply on the variable
          double min() const { return first; }
          /// Lower limit to apply on the variable
          double& min() { return first; }
          /// Upper limit to apply on the variable
          double max() const { return second; }
          /// Upper limit to apply on the variable
          double& max() { return second; }
          double x( double v ) const;
          /// Specify the lower and upper limits on the variable
          void in( double low, double up );
          /// Full variable range allowed
          double range() const;
          /// Have a lower limit?
          bool hasMin() const;
          /// Have an upper limit?
          bool hasMax() const;
          bool passes( double val ) const;
          bool valid() const;

          /// Human-readable expression of the limits
          friend std::ostream& operator<<( std::ostream&, const Limits& );
        private:
          static constexpr double invalid_ = -999.999;
      };
    public:
      Kinematics();
      Kinematics( const Kinematics& kin );
      ~Kinematics();

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
      void dump( std::ostream& os = Logger::get().outputStream ) const;

      /// Incoming particles' momentum (in \f$\text{GeV}/c\f$)
      std::pair<double,double> inp;
      /// Set the incoming particles' momenta (if the collision is symmetric)
      void setSqrtS( double sqrts );
      /// Process centre of mass energy
      double sqrtS() const;
      /// Beam/primary particle's PDG identifier
      std::pair<ParticleCode,ParticleCode> inpdg;
      /// PDG id of the outgoing central particles
      std::vector<ParticleCode> central_system;
      /// Minimum list of central particles required
      std::vector<ParticleCode> minimum_final_state;

      /// Type of kinematics to consider for the phase space
      ProcessMode mode;
      /// Type of structure functions to consider
      StructureFunctions::Type structure_functions;

      struct CutsList {
        CutsList();
        CutsList( const CutsList& cuts );
        /// Cuts on the initial particles kinematics
        std::map<Cuts::InitialState,Limits> initial;
        /// Cuts on the central system produced
        std::map<Cuts::Central,Limits> central;
        std::map<ParticleCode,std::map<Cuts::Central,Limits> > central_particles;
        /// Cuts on the beam remnants system
        std::map<Cuts::Remnants,Limits> remnants;
      };
      CutsList cuts;
  };
}

#endif

