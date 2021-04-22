#ifndef CepGen_Physics_Kinematics_h
#define CepGen_Physics_Kinematics_h

#include "CepGen/Core/ParametersList.h"

#include "CepGen/Physics/Cuts.h"
#include "CepGen/Physics/Modes.h"

#include <iosfwd>
#include <vector>
#include <memory>

namespace cepgen
{
  enum class KTFlux;
  namespace strfun { class Parameterisation; }
  namespace formfac { class Parameterisation; }
  /// List of kinematic constraints to apply on the process phase space.
  class Kinematics
  {
    public:
      Kinematics() = default;
      explicit Kinematics( const ParametersList& );
      ~Kinematics() = default;

      /// Minimal diffractive mass for dissociative proton treatment
      static const double MX_MIN;

      /// List containing all parameters handled
      ParametersList parameters() const;

      /// Set the incoming particles' momenta (if the collision is symmetric)
      Kinematics& setSqrtS( double sqrts );
      /// Process centre of mass energy
      double sqrtS() const;

      /// Set the beams for the type of kinematics to consider
      Kinematics& setMode( const mode::Kinematics& );
      /// Type of kinematics to consider for the phase space
      mode::Kinematics mode() const;

      /// Incoming beams characteristics
      struct Beam
      {
        Beam(); ///< Default constructor
        double pz; ///< Incoming particle momentum, in GeV/c
        pdgid_t pdg; ///< PDG identifier for the beam
        mode::Beam mode; ///< Beam treatment mode
        KTFlux kt_flux; ///< Type of \f$k_{\rm T}\f$-factorised flux to be considered (if any)
      };
      /// Human-readable description of a beam particle/system
      friend std::ostream& operator<<( std::ostream&, const Beam& );

      /// Beam/primary particle's kinematics
      std::pair<Beam,Beam> incoming_beams;
      /// Minimum list of central particles required
      std::vector<pdgid_t> minimum_final_state;

      /// Form factors evaluator
      formfac::Parameterisation* formFactors() const { return form_factors_.get(); }
      /// Set a form factors evaluator object
      Kinematics& setFormFactors( std::unique_ptr<formfac::Parameterisation> );

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
        cuts::Initial initial; ///< Cuts on the initial particles kinematics
        cuts::Central central; ///< Cuts on the central system produced
        PerIdCuts central_particles; ///< Cuts on the central individual particles
        cuts::Remnants remnants; ///< Cuts on the beam remnants system
      } cuts; ///< Phase space cuts
      /// Human-readable description of a full kinematics cuts definition
      friend std::ostream& operator<<( std::ostream&, const CutsList& );

    private:
      /// Type of form factors to consider
      std::shared_ptr<formfac::Parameterisation> form_factors_;
      /// Type of structure functions to consider
      std::shared_ptr<strfun::Parameterisation> str_fun_;
  };
}

#endif
