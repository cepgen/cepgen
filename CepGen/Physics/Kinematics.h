#ifndef CepGen_Physics_Kinematics_h
#define CepGen_Physics_Kinematics_h

#include "CepGen/Core/Hasher.h"

#include "CepGen/Physics/Cuts.h"
#include "CepGen/Physics/HeavyIon.h"

#include <ostream>
#include <vector>
#include <unordered_map>
#include <memory>

namespace CepGen
{
  enum class PDG;
  enum class KTFlux;
  class StructureFunctions;
  /// Type of kinematics to consider for the process
  enum class KinematicsMode
  {
    invalid = -1,
    ElectronProton = 0,     ///< electron-proton elastic case
    ElasticElastic = 1,     ///< proton-proton elastic case
    ElasticInelastic = 2,   ///< proton-proton single-dissociative (or inelastic-elastic) case
    InelasticElastic = 3,   ///< proton-proton single-dissociative (or elastic-inelastic) case
    InelasticInelastic = 4, ///< proton-proton double-dissociative case
    ProtonElectron,
    ElectronElectron
  };
  /// Human-readable format of a process mode (elastic/dissociative parts)
  std::ostream& operator<<( std::ostream&, const KinematicsMode& );
  /// List of kinematic constraints to apply on the process phase space.
  class Kinematics
  {
    public:
      Kinematics();
      ~Kinematics();

      /// Incoming beams characteristics
      struct Beam
      {
        double pz; ///< Incoming particle momentum, in GeV/c
        PDG pdg; ///< PDG identifier for the beam
        KTFlux kt_flux; ///< Type of kT-factorised flux to be considered (if any)
      };
      friend std::ostream& operator<<( std::ostream&, const Beam& );
      /// Beam/primary particle's kinematics
      std::pair<Beam,Beam> incoming_beams;
      /// Set the incoming particles' momenta (if the collision is symmetric)
      void setSqrtS( double sqrts );
      /// Process centre of mass energy
      double sqrtS() const;
      /// Minimum list of central particles required
      std::vector<PDG> minimum_final_state;

      /// Type of kinematics to consider for the phase space
      KinematicsMode mode;
      /// Type of structure functions to consider
      std::shared_ptr<StructureFunctions> structure_functions;

      /// A collection of cuts to apply on the physical phase space
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
      std::string kmr_grid_path;
  };
}

#endif

