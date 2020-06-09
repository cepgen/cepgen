#ifndef CepGen_Physics_Kinematics_h
#define CepGen_Physics_Kinematics_h

#include "CepGen/Core/ParametersList.h"
#include "CepGen/Physics/KinematicsMode.h"
#include "CepGen/Physics/Cuts.h"

#include <iosfwd>
#include <vector>
#include <memory>

namespace cepgen
{
  enum class KTFlux;
  namespace strfun { class Parameterisation; }
  /// List of kinematic constraints to apply on the process phase space.
  class Kinematics
  {
    public:
      Kinematics();
      Kinematics( const ParametersList& );
      ~Kinematics() = default;

      /// Minimal diffractive mass for dissociative proton treatment
      static constexpr double MX_MIN = 1.07; // mp+mpi+-

      /// List containing all parameters handled
      ParametersList parameters() const;

      /// Set the incoming particles' momenta (if the collision is symmetric)
      Kinematics& setSqrtS( double sqrts );
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

      /// Structure functions evaluator
      strfun::Parameterisation* structureFunctions() const { return str_fun_.get(); }
      /// Set a structure functions evaluator object
      Kinematics& setStructureFunctions( std::unique_ptr<strfun::Parameterisation> );
      /// Set the integer-type of structure functions evaluator to build
      Kinematics& setStructureFunctions( int, int );

      /// A collection of cuts to apply on the physical phase space
      struct CutsList
      {
        CutsList();
        InitialCuts initial; ///< Cuts on the initial particles kinematics
        CentralCuts central; ///< Cuts on the central system produced
        PerIdCuts central_particles; ///< Cuts on the central individual particles
        RemnantsCuts remnants; ///< Cuts on the beam remnants system
      } cuts; ///< Phase space cuts
      /// Human-readable description of a full kinematics cuts definition
      friend std::ostream& operator<<( std::ostream&, const CutsList& );

    private:
      /// Type of structure functions to consider
      std::shared_ptr<strfun::Parameterisation> str_fun_;
  };
}

#endif
