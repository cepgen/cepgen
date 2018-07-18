#ifndef CepGen_Physics_Kinematics_h
#define CepGen_Physics_Kinematics_h

#include <iomanip>
#include <algorithm>

#include "CepGen/Core/Logger.h"

#include "CepGen/Physics/ParticleProperties.h"
#include "CepGen/Physics/Cuts.h"
#include "CepGen/Physics/Limits.h"

#include <vector>
#include <unordered_map>

namespace CepGen
{
  enum class PDG;
  class StructureFunctions;
  /// List of kinematic constraints to apply on the process phase space.
  class Kinematics
  {
    public:
      Kinematics();
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
      std::shared_ptr<StructureFunctions> structure_functions;

      struct CutsList
      {
        CutsList();
        template<class T,bool>
        struct hasher
        {
          inline size_t operator()( const T& t ) const {
            return std::hash<T>()( t );
          }
        };
        template<class T>
        struct hasher<T, true>
        {
          inline size_t operator() ( const T& t ) {
            typedef typename std::underlying_type<T>::type enumType;
            return std::hash<enumType>()( static_cast<enumType>( t ) );
          }
        };
        template<class T>
        struct EnumHash
        {
          inline size_t operator()( const T& t ) const {
            return hasher<T,std::is_enum<T>::value>()( t );
          }
        };
        /// Cuts on the initial particles kinematics
        Cuts initial;
        /// Cuts on the central system produced
        Cuts central;
        std::unordered_map<PDG,Cuts,EnumHash<PDG> > central_particles;
        /// Cuts on the beam remnants system
        Cuts remnants;
      };
      CutsList cuts;
  };
}

#endif

