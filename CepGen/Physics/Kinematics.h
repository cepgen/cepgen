#ifndef CepGen_Physics_Kinematics_h
#define CepGen_Physics_Kinematics_h

#include "CepGen/Physics/KinematicsMode.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/Cuts.h"
#include "CepGen/Physics/HeavyIon.h"

#include <ostream>
#include <vector>
#include <unordered_map>
#include <memory>

namespace cepgen
{
  enum class KTFlux;
  class ParametersList;
  namespace strfun { class Parameterisation; }
  /// List of kinematic constraints to apply on the process phase space.
  class Kinematics
  {
    public:
      Kinematics();
      Kinematics( const ParametersList& );
      ~Kinematics() = default;

      /// Set the incoming particles' momenta (if the collision is symmetric)
      void setSqrtS( double sqrts );
      /// Process centre of mass energy
      double sqrtS() const;

      /// Incoming beams characteristics
      struct Beam
      {
        double pz; ///< Incoming particle momentum, in GeV/c
        pdgid_t pdg; ///< PDG identifier for the beam
        KTFlux kt_flux; ///< Type of \f$k_{\rm T}\f$-factorised flux to be considered (if any)
      };
      /// Human-readable description of a beam particle/system
      friend std::ostream& operator<<( std::ostream&, const Beam& );

      /// Beam/primary particle's kinematics
      std::pair<Beam,Beam> incoming_beams;
      /// Minimum list of central particles required
      std::vector<pdgid_t> minimum_final_state;
      /// Type of kinematics to consider for the phase space
      KinematicsMode mode;
      /// Type of structure functions to consider
      std::shared_ptr<strfun::Parameterisation> structure_functions;

      /// A collection of cuts to apply on the physical phase space
      struct CutsList
      {
        CutsList();
        Cuts initial; ///< Cuts on the initial particles kinematics
        Cuts central; ///< Cuts on the central system produced
        PerIdCuts central_particles; ///< Cuts on the central individual particles
        Cuts remnants; ///< Cuts on the beam remnants system
      } cuts; ///< Phase space cuts
  };
}

#endif
