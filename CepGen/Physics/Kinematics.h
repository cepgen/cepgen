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
#include "CepGen/Core/Hasher.h"

namespace CepGen
{
  enum class PDG;
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
      struct HeavyIon {
        static inline HeavyIon Proton() { return HeavyIon{ 1, 1 }; }
        static inline HeavyIon Pb208() { return HeavyIon{ 208, 82 }; }
        static inline HeavyIon fromPDG( const PDG& pdg ) {
          unsigned int ipdg = (unsigned int)pdg;
          return HeavyIon{ (unsigned short)( ipdg % 1000 ), (unsigned short)( ( ipdg / 1000 ) % 1000 ) };
        }
        inline PDG pdg() const { return (PDG)( 1e3*Z+A ); } // (Pythia8 convention/10-1e10)
        /// Mass number
        unsigned short A;
        /// Atomic number
        unsigned short Z;
      };
      friend std::ostream& operator<<( std::ostream&, const HeavyIon& );

      /// Dump all the parameters used in this process cross-section computation
      /// or events generation
      void dump( std::ostream& os = *Logger::get().output ) const;

      struct Beam
      {
        /// Incoming particle's momentum (in \f$\text{GeV}/c\f$)
        double pz;
        PDG pdg;
        HeavyIon hi;
        unsigned short kt_flux;
      };
      /// Beam/primary particle's kinematics
      std::pair<Beam,Beam> incoming_beams;
      /// Set the incoming particles' momenta (if the collision is symmetric)
      void setSqrtS( double sqrts );
      /// Process centre of mass energy
      double sqrtS() const;
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

