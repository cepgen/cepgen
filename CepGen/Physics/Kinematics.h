#ifndef CepGen_Physics_Kinematics_h
#define CepGen_Physics_Kinematics_h

#include <iomanip>
#include <algorithm>

#include "CepGen/Core/Logger.h"
#include "CepGen/StructureFunctions/StructureFunctions.h"

#include "CepGen/Physics/ParticleProperties.h"
#include "CepGen/Physics/Cuts.h"
#include "CepGen/Physics/Limits.h"

#include <vector>
#include <unordered_map>

using std::cout;
using std::string;

namespace CepGen
{
  /// List of kinematic constraints to apply on the process phase space.
  class Kinematics
  {
    public:
      Kinematics();
      Kinematics( const Kinematics& kin );
      Kinematics& operator=( const Kinematics& kin );
      ~Kinematics();

      /// Type of kinematics to consider for the process
      enum class Mode {
        ElectronProton = 0,     ///< electron-proton elastic case
        ElasticElastic = 1,     ///< proton-proton elastic case
        ElasticInelastic = 2,   ///< proton-proton single-dissociative (or inelastic-elastic) case
        InelasticElastic = 3,   ///< proton-proton single-dissociative (or elastic-inelastic) case
        InelasticInelastic = 4, ///< proton-proton double-dissociative case
        ProtonElectron,
        ElectronElectron
      };
      /// Human-readable format of a process mode (elastic/dissociative parts)
      friend std::ostream& operator<<( std::ostream&, const Mode& );

      /// Dump all the parameters used in this process cross-section computation
      /// or events generation
      void dump( std::ostream& os = *Logger::get().output ) const;

      /// Incoming particles' momentum (in \f$\text{GeV}/c\f$)
      std::pair<double,double> inp;
      /// Set the incoming particles' momenta (if the collision is symmetric)
      void setSqrtS( double sqrts );
      /// Process centre of mass energy
      double sqrtS() const;
      /// Beam/primary particle's PDG identifier
      std::pair<PDG,PDG> inpdg;
      /// PDG id of the outgoing central particles
      std::vector<PDG> central_system;
      /// Minimum list of central particles required
      std::vector<PDG> minimum_final_state;

      /// Type of kinematics to consider for the phase space
      Mode mode;
      /// Type of structure functions to consider
      StructureFunctions::Type structure_functions;

      struct CutsList
      {
        CutsList();
        CutsList( const CutsList& cuts );
        CutsList& operator=( const CutsList& cuts );
        struct EnumClassHash
        {
          template <typename T> std::size_t operator()( T t ) const {
            return static_cast<std::size_t>( t );
          }
        };
        /// Cuts on the initial particles kinematics
        /*std::map<Cuts,Limits> initial;
        /// Cuts on the central system produced
        std::map<Cuts,Limits> central;
        std::map<PDG,std::map<Cuts,Limits> > central_particles;
        /// Cuts on the beam remnants system
        std::map<Cuts,Limits> remnants;*/
        /// Cuts on the initial particles kinematics
        std::unordered_map<Cuts,Limits,EnumClassHash> initial;
        /// Cuts on the central system produced
        std::unordered_map<Cuts,Limits,EnumClassHash> central;
        std::map<PDG,std::unordered_map<Cuts,Limits,EnumClassHash> > central_particles;
        //std::unordered_map<PDG,std::unordered_map<Cuts,Limits,EnumClassHash>,PDGHash> central_particles;
        /// Cuts on the beam remnants system
        std::unordered_map<Cuts,Limits,EnumClassHash> remnants;
      };
      CutsList cuts;
  };
}

#endif

